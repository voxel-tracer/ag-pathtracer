#pragma once

#include "camera.h"
#include "material.h"
#include "intersectable.h"

class Scene {
public:
	Scene(vector<shared_ptr<Intersectable>> p, float3 lightpos, float3 lightcolor, const CameraDesc& cam) :
		primitives(move(p)), lightPos(lightpos), lightColor(lightcolor), camera(cam) {}

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
	float3 lightPos;
	float3 lightColor;
	CameraDesc camera;
};

float3 DirectIllumination(const Scene* scene, const float3& I, const float3& N) {
	// check if light source is unobstructed
	Ray ray(I + 0.0001f * N, normalize(scene->lightPos - I));
	Hit hit;
	if (scene->NearestIntersection(ray, hit)) {
		return float3(0); // light source doesn't reach I
	}

	// compute distance to light source
	float r2 = sqrLength(scene->lightPos - I);
	// compute cosine between light direction and normal
	float IoN = clamp(dot(ray.D, N), 0.0f, 1.0f); // notice that D points towards the light
	// account for light intensity
	return scene->lightColor * IoN / r2;
}

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

float3 Trace(const Scene* scene, const Ray& ray, int depth) {
	if (depth == 10) return float3(0.0f);

	Hit hit;
	if (!scene->NearestIntersection(ray, hit)) {
		// Background color
		// RGB 173, 216, 229
		// 0-1 .68, .85, .89
		//return rgb2lin({ .68f, .85f, 0.89f });
		return { .4f, .45f, .5f };
	}
	const Material* mat = hit.mat;

	if (mat->type == DIFFUSE) {
		// TODO material should just evaluate the color
		return mat->baseColor->value(hit.u, hit.v) * DirectIllumination(scene, hit.I, hit.N);
	}
	if (mat->type == MIRROR) {
		// TODO material should compute reflected direction
		float3 D = reflect(ray.D, hit.N);
		float3 O = hit.I + 0.0001f * D;
		return mat->baseColor->value(hit.u, hit.v) * Trace(scene, Ray(O, D), depth + 1);
	}
	if (mat->type == GLASS) {
		// figure out if we are entering or exiting the object
		bool entering = dot(ray.D, hit.N) < 0;
		float3 N = entering ? hit.N : -hit.N; // make sure surface normal is on same side as ray

		// absorption can only happen inside a glass (entering == false)
		float3 transmission(1.0f);
		float3 matColor = mat->baseColor->value(hit.u, hit.v);
		if (!entering && length(matColor) > 0.0f) {
			//transmission.x = exp(-hit.color.x * ray.t);
			//transmission.y = exp(-hit.color.y * ray.t);
			//transmission.z = exp(-hit.color.z * ray.t);
			transmission.x = pow(matColor.x, ray.t);
			transmission.y = pow(matColor.y, ray.t);
			transmission.z = pow(matColor.z, ray.t);
		}

		// use it to compute refraction ratio assuming external material is air
		float etaI, etaT;
		if (entering) { // air -> material
			etaI = 1.0f; // air
			etaT = mat->ref_idx;
		}
		else { // material -> air
			etaI = mat->ref_idx;
			etaT = 1.0f;
		}

		// thetaI is the incident angle
		// notice that we are using ray.D and the normal facing N so cos_thetaI will always be negative
		float cos_thetaI = fmaxf(dot(ray.D, N), -1.0f);

		// compute how much light is reflected
		float Fr = 1.0f;
		float3 T;
		float3 color;
		if (refract(ray.D, N, cos_thetaI, etaI / etaT, T)) {
			// some light is refracted, compute Fresnel reflection
			Fr = FresnelSchlick(etaI, etaT, cos_thetaI);

			// trace refracted ray
			color = (1 - Fr) * transmission * Trace(scene, Ray(hit.I + .0001f * (-N), T), depth + 1);
		}

		// trace reflected ray
		float3 R = reflect(ray.D, N);
		color += Fr * transmission * Trace(scene, Ray(hit.I + .0001f * N, R), depth + 1);

		return color;
	}

	return make_float3(0.0f);
}
