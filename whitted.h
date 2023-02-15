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
	virtual float3 Li(const Ray& ray, const Scene& scene, int depth = 0) const = 0;

protected:
	float3 SpecularReflection(const Ray& ray, const Hit& hit, const Scene& scene, int depth) const {
		if (isblack(hit.specular)) return float3(0.f);
			
		float3 D = reflect(ray.D, hit.N);
		float3 O = hit.I + EPSILON * D;
		return hit.specular * Li(Ray(O, D), scene, depth + 1);
	}

	float3 SpecularTransmission(const Ray& ray, const Hit& hit, const Scene& scene, int depth) const {
		if (isblack(hit.transmission)) return float3(0.f);

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
			L += (1 - Fr) * transmission * Li(Ray(hit.I + EPSILON * (-N), T), scene, depth + 1);
		}

		// TODO we should also account for reflections (total internal reflection + Fresnel reflection at the surface)
		// trace reflected ray
		// float3 R = reflect(ray.D, N);
		// L += Fr * transmission * Li(Ray(hit.I + EPSILON * N, R), scene, depth + 1);

		return L;
	}
};

class WhittedIntegrator : public Integrator {
public:
	WhittedIntegrator(int maxDepth = 5) : maxDepth(maxDepth) {}

	virtual float3 Li(const Ray& ray, const Scene& scene, int depth = 0) const override {
		float3 L(0.f);

		Hit hit;
		if (!scene.NearestIntersection(ray, hit)) {
			for (const auto& light : scene.lights) L += light->Le(ray);
			return L;
		}

		hit.EvalMaterial();

		if (!isblack(hit.diffuse)) 
			L += hit.diffuse * DirectIllumination(scene, hit.I, hit.N);

		if (depth + 1 < maxDepth) {
			L += SpecularReflection(ray, hit, scene, depth);
			L += SpecularTransmission(ray, hit, scene, depth);
		}

		return L;
	}
protected:
	float3 DirectIllumination(const Scene& scene, const float3 I, const float3& N) const {
		float3 L(0.f);
		float3 wi;
		Ray shadowRay;
		Hit tmpHit;

		for (const auto& light : scene.lights) {
			auto Li = light->Sample_Li(I, &wi, &shadowRay);
			if (isblack(Li)) continue; // light didn't return any radiance toward I

			auto WiDotN = dot(wi, N);
			if (WiDotN <= 0.f) continue; // light is behind the surface

			// TODO we need an IntersectP() that just returns true if any intersection is found
			if (scene.NearestIntersection(shadowRay, tmpHit)) continue; // light is obstructed

			L += Li * WiDotN;
		}

		return L;
	}

	int maxDepth;
};