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

		for (const auto& light : scene.lights) {
			float3 S, lightN;
			light->Sample(I, &S, &lightN);
			float3 Wi = S - I;
			float dist = length(Wi);
			Wi /= dist;

			float cos_o = dot(-Wi, lightN);
			float cos_i = dot(Wi, N);
			if (cos_o <= 0 || cos_i <= 0.f) continue; // light is behind the surface
			// trace a shadow ray
			Ray shadowRay(I + EPSILON * Wi, Wi, dist - EPSILON);
			Hit tmpHit;
			if (scene.NearestIntersection(shadowRay, tmpHit)) continue; // light is obstructed
			float3 Li = light->L / (dist * dist);
			L += Li * cos_i;
		}

		return L;
	}

	int maxDepth;
};

class SimplePT : public Integrator {
public:
	SimplePT(int maxDepth = 5) : MaxDepth(maxDepth) {}

	virtual float3 Li(const Ray& ray, const Scene& scene, int depth = 0) const override {
		float3 L(0.f);

		Hit hit;
		if (!scene.NearestIntersection(ray, hit)) {
			// Ray left the scene, handle infinite lights
			for (const auto& light : scene.lights) L += light->Le(ray);
			return L;
		}

		hit.EvalMaterial();

		// terminate if we hit a light source
		if (!isblack(hit.emission)) return hit.emission;
		// terminate if we exceed MaxDepth
		if (depth + 1 > MaxDepth) return L;
		
		// continue in random direction
		auto R = DiffuseReflection(hit.N);
		Ray newRay(hit.I + EPSILON * R, R);
		// update throughput
		auto BRDF = hit.diffuse * INVPI; // diffuse brdf = albedo / pi
		auto Ei = Li(newRay, scene, depth + 1) * dot(hit.N, R); // irradiance
		return TWOPI * BRDF * Ei;
	}

protected:
	float3 DiffuseReflection(const float3& N) const {
		return RandomInHemisphere(N);
	}

	int MaxDepth;
};


class SimplePT2 : public Integrator {
public:
	SimplePT2(int maxDepth = 5) : MaxDepth(maxDepth) {}

	virtual float3 Li(const Ray& ray, const Scene& scene, int depth = 0) const override {
		float3 L(0.f);

		Hit hit;
		if (!scene.NearestIntersection(ray, hit)) {
			// Ray left the scene, handle infinite lights
			for (const auto& light : scene.lights) L += light->Le(ray);
			return L;
		}

		hit.EvalMaterial();

		// terminate if we hit a light source
		if (!isblack(hit.emission)) return float3(0.f);

		// pick one random light
		int lights = (int)(scene.lights.size());
		int lightIdx = clamp((int)(RandomFloat() * lights), 0, lights - 1);
		const auto& light = scene.lights[lightIdx];

		// sample point on the light
		float3 S, lightN;
		if (!light->Sample(hit.I, &S, &lightN)) return false; // light doesn't support sampling
		float3 toL = S - hit.I;
		float dist = length(toL);
		toL /= dist;
		float cos_o = dot(-toL, lightN);
		float cos_i = dot(toL, hit.N);
		if (cos_i <= 0 || cos_o <= 0) return float3(0.f);
		// light is not behind surface point, trace shadow ray
		Ray newRay(hit.I + EPSILON * toL, toL, dist - 2 * EPSILON);
		if (scene.NearestIntersection(newRay, hit)) return float3(0.f); // occluded light
		// light is visible, calculate transport
		auto BRDF = hit.diffuse * INVPI; // diffuse brdf = albedo / pi
		float solidAngle = (cos_o * light->Area()) / (dist * dist);
		return BRDF * lights * light->L * solidAngle * cos_i;
	}

protected:
	float3 DiffuseReflection(const float3& N) const {
		return RandomInHemisphere(N);
	}

	int MaxDepth;
};