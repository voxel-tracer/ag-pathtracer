#pragma once

#include "camera.h"
#include "material.h"
#include "intersectable.h"
#include "lights.h"

class Scene {
public:
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
	virtual float3 Li(const Ray& ray, const Scene& scene, int depth = 0) const = 0;
};

inline float PowerHeuristic(int nf, float fPdf, int ng, float gPdf) {
	float f = nf * fPdf, g = ng * gPdf;
	return (f * f) / (f * f + g * g);
}

inline float3 EstimateDirect(const SurfaceInteraction& si, const float2& uScattering, const Light& light, const float2& uLight, const Scene& scene) {
	float3 Ld(0.f);
	float3 wi;
	float lightPdf = 0, scatteringPdf = 0;
	VisibilityTester visibility;
	float3 Li = light.Sample_Li(si, uLight, &wi, &lightPdf, &visibility);
	if (lightPdf > 0 && !IsBlack(Li)) {
		// Evaluate BSDF for light sampling strategy
		float3 f = si.bsdf.f(si.wo, wi, true) * absdot(wi, si.shading.n);
		scatteringPdf = si.bsdf.Pdf(si.wo, wi, true);
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
	{
		float3 f;
		// Sample scattered direction
		f = si.bsdf.Sample_f(si.wo, &wi, uScattering, &scatteringPdf, true);
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

inline float3 UniformSampleOneLight(const SurfaceInteraction& si, const Scene& scene) {
	// Randomly choose a single light to sample, _light_
	int nLights = int(scene.lights.size());
	if (nLights == 0) return float3(0.f);
	int numLight = std::min((int)(RandomFloat() * nLights), nLights - 1);
	float lightPdf = 1.f / nLights;
	const std::shared_ptr<Light>& light = scene.lights[numLight];
	float2 uLight(RandomFloat(), RandomFloat());
	float2 uScattering(RandomFloat(), RandomFloat());
	return EstimateDirect(si, uScattering, *light, uLight, scene) / lightPdf;
}

class DbgIntegrator : public Integrator {
public:
	float3 Li(const Ray& ray, const Scene& scene, int depth = 0) const {
		SurfaceInteraction si;
		if (scene.NearestIntersection(ray, si)) {
			if (si.uv.x == 0 || si.uv.y == 0)
				return float3(1, 0, 0);
			return (float3(si.uv.x, si.uv.y, 0)) / 5;
		}
		return float3(0.f);
	}
};

class PathTracer : public Integrator {
public:
	PathTracer(int maxDepth = 5) : MaxDepth(maxDepth) {}

	float3 Li(const Ray& r, const Scene& scene, int depth = 0) const {
		float3 beta(1.f); // current ray throughput
		float3 L(0.f);

		Ray ray = r;
		bool specularBounce = false;
		int bounces;

		for (bounces = 0;; bounces++) {

			// Intersect _ray_ with scene and store intersection in _isect_
			SurfaceInteraction isect;
			bool foundIntersection = scene.NearestIntersection(ray, isect);

			// Possibly add emitted light at intersection
			if (bounces == 0 || specularBounce) {
				// Add emitted light at path vertex or from the environment
				if (foundIntersection) {
					L += beta * isect.Le(-ray.D);
				} else {
					for (const auto& light : scene.lights)
						L += beta * light->Le(ray);
				}
			}

			// Terminate ray if ray escaped or _maxDepth_ reached
			if (!foundIntersection || bounces >= MaxDepth) break;

			if (!isect.EvalMaterial()) {
				// Skip intersection due to null BSDF (most likely we hit a light source)
				// It basically ignores the intersection all together
				
				// IMPORTANT: because I don't handle single sided area lights, when a specular ray hits a sphere it may collect energy from both
				// sides of the sphere which is not something we want
				ray = Ray(isect.p + EPSILON * ray.D, ray.D);
				bounces--;
				continue;
			}

			// Sample illumination from lights to find path contribution.
			// (But skip for perfectly specular BSDFs).
			if (!isect.bsdf.IsPerfectlySpecular()) {
				L += beta * UniformSampleOneLight(isect, scene);
			}

			// Sample BSDF to get new path direction
			float3 wo = -ray.D, wi;
			float2 u(RandomFloat(), RandomFloat());
			float pdf;
			bool sampledSpecular = false;
			float3 f = isect.bsdf.Sample_f(wo, &wi, u, &pdf, false, &sampledSpecular);
			if (IsBlack(f) || pdf == 0) break;
			beta *= f * absdot(wi, isect.shading.n) / pdf;
			specularBounce = sampledSpecular;

			// Russian roulette
			float maxComponent = max(beta.x, max(beta.y, beta.z));
			if (maxComponent < 1 && depth > 3) {
				float q = std::max(.05f, 1 - maxComponent);
				if (RandomFloat() < q) break;
				beta /= 1 - q;
			}

			ray = Ray(isect.p + EPSILON * wi, wi);
		}

		return L;
	}

protected:

	int MaxDepth;
};