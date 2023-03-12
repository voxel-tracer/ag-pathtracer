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

bool refract(const float3& I, const float3& N, float cos_thetaI, float etai_over_etat, float3& T) {
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

float FresnelSchlick(float etaI, float etaT, float cos_thetaI) {
	float r0 = (etaI - etaI) / (etaI + etaI);
	r0 = r0 * r0;
	return r0 + (1 - r0) * powf((1 + cos_thetaI), 5); // we use +cos_thetaI instead of -cos_thetaI because we used ray.D to compute the cosine
}

class Integrator {
public:
	virtual float3 Li(const Ray& ray, const Scene& scene, int depth = 0, bool isSpecular = false) const = 0;
};

class PathTracer : public Integrator {
public:
	PathTracer(int maxDepth = 5) : MaxDepth(maxDepth) {}

	virtual float3 Li(const Ray& ray, const Scene& scene, int depth = 0, bool isSpecular = false) const override {
		float3 T(1.f); // current ray throughput
		float3 E(0.f);

		Ray curRay = ray;
		while (true) {
			SurfaceInteraction hit;
			if (!scene.NearestIntersection(curRay, hit)) {
				// Ray left the scene, handle infinite lights
				for (const auto& light : scene.lights) E += T * light->Le(curRay);
				break;
			}

			hit.EvalMaterial();

			// terminate if we hit a light source
			if (!isblack(hit.emission)) {
				if (isSpecular) E += T * hit.emission;
				break;
			}

			E += T * EstimateDirect(scene, hit);

			// terminate if we exceed MaxDepth
			if (depth + 1 > MaxDepth) break;

			// for now assume a material can only be diffuse or specular or refractive
			float3 R;
			if (hit.hasBSDF) {
				T *= SampleMicrofacet(hit, scene, &R);
				isSpecular = false;
			}
			else if (!isblack(hit.diffuse)) {
				T *= SampleIndirectLight(hit, scene, &R);
				isSpecular = false;
			}

			// Russian roulette
			float p = clamp(max(T.x, max(T.y, T.z)), EPSILON, 1.f);
			if (RandomFloat() > p) break;
			T *= 1.f / p; // add the energy we lose by randomly killing paths

			curRay = Ray(hit.p + EPSILON * R, R);
			depth++;
		}

		return E;
	}

protected:

	float3 SampleMicrofacet(const SurfaceInteraction& hit, const Scene& scene, float3* R) const {
		float pdf;
		float2 u(RandomFloat(), RandomFloat());
		float3 BRDF = hit.bsdf.Sample_f(hit.wo, R, u, &pdf);
		if (pdf == 0) return float3(0.f);
		return BRDF * absdot(hit.shading.n, *R) / pdf;
	}

	float3 SampleDirectionInHemisphere(const float3& N, float* pdf) const {
		auto R = CosineWeightedRandomInHemisphere(N);
		*pdf = dot(N, R) * INVPI;
		return R;
	}

	float3 SampleIndirectLight(const SurfaceInteraction& hit, const Scene& scene, float3* R) const {
		float pdf;
		*R = SampleDirectionInHemisphere(hit.n, &pdf);
		// update throughput
		auto BRDF = hit.diffuse * INVPI; // diffuse brdf = albedo / pi
		return BRDF * dot(hit.n, *R) / pdf;
	}

	float3 IndirectLight(const SurfaceInteraction& hit, const Scene& scene, int depth) const {
		float pdf;
		auto R = SampleDirectionInHemisphere(hit.n, &pdf);
		Ray newRay(hit.p + EPSILON * R, R);
		// update throughput
		auto BRDF = hit.diffuse * INVPI; // diffuse brdf = albedo / pi
		auto Ei = Li(newRay, scene, depth + 1) * dot(hit.n, R); // irradiance
		return BRDF * Ei / pdf;
	}

	float3 EstimateDirect(const Scene& scene, const SurfaceInteraction& hit) const {
		if (isblack(hit.diffuse) && !hit.hasBSDF) return float3(0.f);

		// pick one random light
		int lights = (int)(scene.lights.size());
		int lightIdx = clamp((int)(RandomFloat() * lights), 0, lights - 1);
		const auto& light = scene.lights[lightIdx];
		float3 S, lightN;
		if (!light->Sample(hit.p, &S, &lightN)) return false; // light doesn't support sampling
		float3 toL = S - hit.p;
		float dist = length(toL);
		toL /= dist;
		float cos_o = dot(-toL, lightN);
		float cos_i = absdot(toL, hit.shading.n);
		if (cos_i <= 0 || cos_o <= 0) return float3(0.f);
		// light is not behind surface point, trace shadow ray
		Ray newRay(hit.p + EPSILON * toL, toL, dist - 2 * EPSILON);
		SurfaceInteraction tmpHit;
		if (scene.NearestIntersection(newRay, tmpHit)) return float3(0.f); // occluded light
		// light is visible, calculate transport
		float3 BRDF;
		if (hit.hasBSDF) {
			BRDF = hit.bsdf.f(hit.wo, toL);
		}
		else {
			BRDF = hit.diffuse * INVPI; // diffuse brdf = albedo / pi
		}
		float solidAngle = (cos_o * light->Area()) / (dist * dist);
		return BRDF * (float)lights * light->L * solidAngle * cos_i;
	}

	int MaxDepth;
};