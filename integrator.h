#pragma once

#include "camera.h"
#include "material.h"
#include "intersectable.h"
#include "lights.h"

class Scene {
public:
	Scene(vector<shared_ptr<Intersectable>> p, const CameraDesc& cam) :
		primitives(move(p)), camera(cam) {}

	bool NearestIntersection(const Ray& ray, SurfaceInteraction& hit) const {
		bool found = false;
		for (auto& primitive : primitives) {
			if (primitive->Intersect(ray, hit)) {
				found = true;
			}
		}
		return found;
	}

	vector<shared_ptr<Intersectable>> primitives;
	vector<shared_ptr<Light>> lights;
	CameraDesc camera;
};

inline bool refract(const float3& I, const float3& N, float cos_thetaI, float etai_over_etat, float3& T) {
	// thetaI is the incident angle
	// N is always facing the ray origin and I is pointing toward the surface => cos_thetaI will always be negative
	float sin_thetaI = sqrtf(1.f - cos_thetaI * cos_thetaI);
	if (sin_thetaI > (1 / etai_over_etat))
		return false; // total internal reflection

	// compute refracted direction
	float sin_thetaT2 = etai_over_etat * etai_over_etat * sin_thetaI * sin_thetaI;
	T = etai_over_etat * I - (etai_over_etat * cos_thetaI + sqrtf(1 - sin_thetaT2)) * N;
	return true;
}

inline float FresnelSchlick(float etaI, float etaT, float cos_thetaI) {
	float r0 = (etaI - etaI) / (etaI + etaI);
	r0 = r0 * r0;
	return r0 + (1 - r0) * powf((1 + cos_thetaI), 5); // we use +cos_thetaI instead of -cos_thetaI because we used ray.D to compute the cosine
}

class Integrator {
public:
	virtual float3 Li(const Ray& ray, const Scene& scene, int depth = 0, bool isSpecular = true) const = 0;
};

inline float PowerHeuristic(int nf, float fPdf, int ng, float gPdf) {
	float f = nf * fPdf, g = ng * gPdf;
	return (f * f) / (f * f + g * g);
}

inline float3 EstimateDirect(const SurfaceInteraction& si, const float2& uScattering, const Light& light, const float2& uLight, const Scene& scene, bool MIS = false) {
	float3 Ld(0.f);
	float3 wi;
	float lightPdf = 0, scatteringPdf = 0;
	VisibilityTester visibility;
	float3 Li = light.Sample_Li(si, uLight, &wi, &lightPdf, &visibility);
	if (lightPdf > 0 && !IsBlack(Li)) {
		// Evaluate BSDF for light sampling strategy
		float3 f = si.bsdf.f(si.wo, wi) * absdot(wi, si.shading.n);
		scatteringPdf = si.bsdf.Pdf(si.wo, wi);
		if (!IsBlack(f)) {
			// Compute effect of visibility for light source sample
			if (!visibility.Unoccluded(scene))
				Li = float3(0.f);

			// Add light's contribution to reflected radiance
			if (!IsBlack(Li)) {
				// TODO handle delta lights (point light)
				float weight = PowerHeuristic(1, lightPdf, 1, scatteringPdf);
				Ld += f * Li * weight / lightPdf;
			}
		}
	}

	// Sample BSDF with multiple importance sampling
	if (MIS)
	{
		float3 f;
		// Sample scattered direction
		f = si.bsdf.Sample_f(si.wo, &wi, uScattering, &scatteringPdf);
		f *= absdot(wi, si.shading.n);

		if (!isblack(f) && scatteringPdf > 0) {
			// Account for light contribution along sampled direction _wi_
			float weight = 1;
			lightPdf = light.Pdf_Li(si, wi);
			if (lightPdf == 0) return Ld;
			weight = PowerHeuristic(1, scatteringPdf, 1, lightPdf);

			// Find Intersection
			SurfaceInteraction lightIsect;
			Ray ray(si.p + EPSILON * wi, wi);
			bool foundSurfaceInteraction = scene.NearestIntersection(ray, lightIsect);
			// Add light contribution from material sampling
			float3 Li(0.f);
			if (foundSurfaceInteraction) {
				if (lightIsect.shape->GetAreaLight() == &light)
					Li = light.L; // TODO we should use Light::L(si, -wi)
			}
			else
				Li = light.Le(ray);
			if (!isblack(Li)) Ld += f * Li * weight / scatteringPdf;
		}
	}

	return Ld;
}

inline float3 UniformSampleOneLight(const SurfaceInteraction& si, const Scene& scene, bool MIS = false) {
	// Randomly choose a single light to sample, _light_
	int nLights = int(scene.lights.size());
	if (nLights == 0) return float3(0.f);
	int numLight = std::min((int)(RandomFloat() * nLights), nLights - 1);
	float lightPdf = 1.f / nLights;
	const std::shared_ptr<Light>& light = scene.lights[numLight];
	float2 uLight(RandomFloat(), RandomFloat());
	float2 uScattering(RandomFloat(), RandomFloat());
	return EstimateDirect(si, uScattering, *light, uLight, scene, MIS) / lightPdf;
}

class PathTracer : public Integrator {
public:
	PathTracer(int maxDepth = 5, bool MIS = true) : MaxDepth(maxDepth), MIS(MIS) {}

	virtual float3 Li(const Ray& ray, const Scene& scene, int depth = 0, bool isSpecular = true) const override {
		float3 T(1.f); // current ray throughput
		float3 E(0.f);

		Ray curRay = ray;
		while (true) {
			SurfaceInteraction hit;
			if (!scene.NearestIntersection(curRay, hit)) {
				// Ray left the scene, handle infinite lights
				if (isSpecular)
					for (const auto& light : scene.lights) 
						E += T * light->Le(curRay);
				break;
			}

			hit.EvalMaterial();

			// terminate if we hit a light source
			if (!isblack(hit.emission)) {
				if (isSpecular) E += T * hit.emission;
				break;
			}

			E += T * UniformSampleOneLight(hit, scene, MIS);

			// terminate if we exceed MaxDepth
			if (depth + 1 > MaxDepth) break;

			// for now assume a material can only be diffuse or specular or refractive
			float3 wi;
			if (hit.hasBSDF) {
					T *= SampleMicrofacet(hit, scene, &wi);
				isSpecular = false;
			}

			// Russian roulette
			float p = clamp(max(T.x, max(T.y, T.z)), EPSILON, 1.f);
			if (RandomFloat() > p) break;
			T *= 1.f / p; // add the energy we lose by randomly killing paths

			curRay = Ray(hit.p + EPSILON * wi, wi);
			depth++;
		}

		return E;
	}

protected:

	float3 SampleMicrofacet(const SurfaceInteraction& hit, const Scene& scene, float3* wi) const {
		float pdf;
		float2 u(RandomFloat(), RandomFloat());
		float3 BRDF = hit.bsdf.Sample_f(hit.wo, wi, u, &pdf);
		if (IsBlack(BRDF) || pdf == 0) return float3(0.f);
		return BRDF * absdot(hit.shading.n, *wi) / pdf;
	}

	int MaxDepth;
	bool MIS;
};