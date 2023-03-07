#pragma once

#include "camera.h"
#include "material.h"
#include "intersectable.h"
#include "lights.h"

class Scene {
public:
	Scene(vector<shared_ptr<Intersectable>> p, const CameraDesc& cam) :
		primitives(move(p)), camera(cam) {}

	bool NearestIntersection(const Ray& ray, Hit& hit) const {
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

protected:
	float3 SpecularReflection(const Ray& ray, const Hit& hit, const Scene& scene, int depth) const {
		float3 D = reflect(ray.D, hit.N);
		float3 O = hit.I + EPSILON * D;
		return hit.specular * Li(Ray(O, D), scene, depth + 1, true);
	}

	float3 SpecularTransmission(const Ray& ray, const Hit& hit, const Scene& scene, int depth) const {
		// figure out if we are entering or exiting the object
		auto entering = dot(ray.D, hit.N) < 0;
		auto N = entering ? hit.N : -hit.N; // make sure surface normal is on same side as ray

		// absorption can only happen inside a glass (entering == false)
		float3 transmission(1.0f);
		if (!entering) transmission = pow(hit.transmission, ray.t);

		// compute refraction index assuming external material is air
		float etaI, etaT;
		if (entering) { // air -> material
			etaI = 1.0f; // air
			etaT = hit.mat->ref_idx;
		}
		else { // material -> air
			etaI = hit.mat->ref_idx;
			etaT = 1.0f;
		}

		// thetaI is the incident angle
		// notice that we are using ray.D and the normal facing N so cos_thetaI will always be negative
		float cos_thetaI = fmaxf(dot(ray.D, N), -1.0f);

		float3 T;
		float Fr = 1;
		float3 L(0.f);
		if (refract(ray.D, N, cos_thetaI, etaI / etaT, T)) {
			// some light is refracted, compute Fresnel reflection
			Fr = FresnelSchlick(etaI, etaT, cos_thetaI);
			// trace refracted ray
			L += (1 - Fr) * transmission * Li(Ray(hit.I + EPSILON * (-N), T), scene, depth + 1, true);
		}

		// TODO we should also account for reflections (total internal reflection + Fresnel reflection at the surface)
		// trace reflected ray
		// float3 R = reflect(ray.D, N);
		// L += Fr * transmission * Li(Ray(hit.I + EPSILON * N, R), scene, depth + 1);

		return L;
	}
};

class PathTracer : public Integrator {
public:
	PathTracer(int maxDepth = 5) : MaxDepth(maxDepth) {}

	virtual float3 Li(const Ray& ray, const Scene& scene, int depth = 0, bool isSpecular = false) const override {
		float3 L(0.f);

		Hit hit;
		if (!scene.NearestIntersection(ray, hit)) {
			// Ray left the scene, handle infinite lights
			for (const auto& light : scene.lights) L += light->Le(ray);
			return L;
		}

		hit.EvalMaterial();

		// terminate if we hit a light source
		if (!isblack(hit.emission)) 
			return isSpecular ? hit.emission : float3(0.f);

		L += DirectLight(scene, -ray.D, hit);
		// terminate if we exceed MaxDepth
		if (depth + 1 > MaxDepth) return L;

		if (depth + 1 < MaxDepth) {
			// for now assume a material can only be diffuse or specular or refractive
			if (!isblack(hit.diffuse))
				L += IndirectLight(hit, scene, depth);
			else if (!isblack(hit.transmission))
				L += SpecularTransmission(ray, hit, scene, depth);
			else if (!isblack(hit.specular))
				L += SpecularReflection(ray, hit, scene, depth);
		}

		return L;
	}

protected:
	virtual float3 SampleDirectionInHemisphere(const float3& N, float* pdf) const {
		auto R = CosineWeightedRandomInHemisphere(N);
		*pdf = dot(N, R) * INVPI;
		return R;
	}

	float3 IndirectLight(const Hit& hit, const Scene& scene, int depth) const {
		float pdf;
		auto R = SampleDirectionInHemisphere(hit.N, &pdf);
		Ray newRay(hit.I + EPSILON * R, R);
		// update throughput
		auto BRDF = hit.diffuse * INVPI; // diffuse brdf = albedo / pi
		auto Ei = Li(newRay, scene, depth + 1) * dot(hit.N, R); // irradiance
		return BRDF * Ei / pdf;
	}

	float3 DirectLight(const Scene& scene, const float3& wo, const Hit& hit) const {
		if (isblack(hit.diffuse) && !hit.hasBSDF) return float3(0.f);

		// pick one random light
		int lights = (int)(scene.lights.size());
		int lightIdx = clamp((int)(RandomFloat() * lights), 0, lights - 1);
		const auto& light = scene.lights[lightIdx];
		float3 S, lightN;
		if (!light->Sample(hit.I, &S, &lightN)) return false; // light doesn't support sampling
		float3 toL = S - hit.I;
		float dist = length(toL);
		toL /= dist;
		float cos_o = dot(-toL, lightN);
		float cos_i = absdot(toL, hit.ShadingN);
		if (cos_i <= 0 || cos_o <= 0) return float3(0.f);
		// light is not behind surface point, trace shadow ray
		Ray newRay(hit.I + EPSILON * toL, toL, dist - 2 * EPSILON);
		Hit tmpHit;
		if (scene.NearestIntersection(newRay, tmpHit)) return float3(0.f); // occluded light
		// light is visible, calculate transport
		float3 BRDF;
		if (hit.hasBSDF) {
			BRDF = hit.bsdf.f(wo, toL);
		}
		else {
			BRDF = hit.diffuse * INVPI; // diffuse brdf = albedo / pi
		}
		float solidAngle = (cos_o * light->Area()) / (dist * dist);
		return BRDF * (float)lights * light->L * solidAngle * cos_i;
	}

	int MaxDepth;
};

class PathTracer2 : public PathTracer {
public:
	PathTracer2(int maxDepth = 5) : PathTracer(maxDepth) {}


	virtual float3 Li(const Ray& ray, const Scene& scene, int depth = 0, bool isSpecular = false) const override {
		float3 T(1.f); // current ray throughput
		float3 E(0.f);

		Ray curRay = ray;
		while (true) {
			Hit hit;
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

			float3 wo = -curRay.D;

			E += T * DirectLight(scene, wo, hit);

			// terminate if we exceed MaxDepth
			if (depth + 1 > MaxDepth) break;

			// for now assume a material can only be diffuse or specular or refractive
			float3 R;
			if (hit.hasBSDF) {
				T *= SampleMicrofacet(hit, scene, wo, &R);
				isSpecular = false;
			} else if (!isblack(hit.diffuse)) {
				T *= SampleIndirectLight(hit, scene, &R);
				isSpecular = false;
			}
			else if (!isblack(hit.specular) || !isblack(hit.transmission)) {
				auto entering = dot(curRay.D, hit.N) < 0;
				auto N = entering ? hit.N : -hit.N; // make sure surface normal is on same side as ray

				// absorption can only happen inside a glass (entering == false)
				float3 transmission(1.0f);
				if (!entering) transmission = pow(hit.transmission, curRay.t);

				// compute refraction index assuming external material is air
				float etaI, etaT;
				if (entering) { // air -> material
					etaI = 1.0f; // air
					etaT = hit.mat->ref_idx;
				}
				else { // material -> air
					etaI = hit.mat->ref_idx;
					etaT = 1.0f;
				}

				// thetaI is the incident angle
				// notice that we are using ray.D and the normal facing N so cos_thetaI will always be negative
				float cos_thetaI = fmaxf(dot(curRay.D, N), -1.0f);

				// Compute Fresnel reflection
				auto Fr = FresnelSchlick(etaI, etaT, cos_thetaI);

				// pick reflection or refraction randomly according to the Fresnel reflection
				if (RandomFloat() < Fr) {
					// Specular reflection
					R = reflect(curRay.D, N);
					T *= hit.specular;
					isSpecular = true;
				}
				else {
					// Specular refraction
					if (!refract(curRay.D, N, cos_thetaI, etaI / etaT, R)) {
						// total internal reflection
						break;
					}
					
					T *= transmission;
					isSpecular = true;
				}
			}

			// Russian roulette
			float p = clamp(max(T.x, max(T.y, T.z)), EPSILON, 1.f);
			if (RandomFloat() > p) break;
			T *= 1.f / p; // add the energy we lose by randomly killing paths

			curRay = Ray(hit.I + EPSILON * R, R);
			depth++;
		}

		return E;
	}

protected:
	float3 SampleIndirectLight(const Hit& hit, const Scene& scene, float3* R) const {
		float pdf;
		*R = SampleDirectionInHemisphere(hit.N, &pdf);
		// update throughput
		auto BRDF = hit.diffuse * INVPI; // diffuse brdf = albedo / pi
		return BRDF * dot(hit.N, *R) / pdf;
	}

	float3 SampleMicrofacet(const Hit& hit, const Scene& scene, const float3& wo, float3* R) const {
		float pdf;
		float2 u(RandomFloat(), RandomFloat());
		float3 BRDF = hit.bsdf.Sample_f(wo, R, u, &pdf);
		if (pdf == 0) return float3(0.f);
		return BRDF * absdot(hit.ShadingN, *R) / pdf;
	}
};