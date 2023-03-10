#pragma once

#include "material.h"
#include "camera.h"
#include "reflection.h"

class BSDF;

class Hit {
public:
	float3 I; // intersection point
	float3 N; // geometric normal at intersection
	float3 ShadingN; // shading normal at intersection

	const Material* mat;
	float u;
	float v;

	float3 dpdu = float3(0.f);
	float3 dpdv = float3(0.f);

	float3 diffuse;
	float3 specular;
	float3 transmission;
	float3 emission;
	BSDF bsdf;
	bool hasBSDF = false;

	void SetShadingGeometry(const float3& dpdus, const float3& dpdvs, bool orientationIsAuthorative) {
		// Compute shading normal for surface interaction
		ShadingN = normalize(cross(dpdus, dpdvs));
		if (orientationIsAuthorative)
			N = Faceforward(N, ShadingN);
		else
			ShadingN = Faceforward(ShadingN, N);
		dpdu = dpdus;
		dpdv = dpdvs;
	}

	void EvalMaterial() {
		if (mat->microfacet) {
			bsdf = BSDF(N, ShadingN, dpdu);
			bsdf.AddBxDF(mat->microfacet.get());
			hasBSDF = true;
		}
		else if (mat->disney) {
			bsdf = BSDF(N, ShadingN, dpdu);
			if (mat->disney->diffuse)
				bsdf.AddBxDF(mat->disney->diffuse.get());
			if (mat->disney->retro)
				bsdf.AddBxDF(mat->disney->retro.get());
			if (mat->disney->microfacet)
				bsdf.AddBxDF(mat->disney->microfacet.get());
			hasBSDF = true;
		}

		diffuse = (mat->diffuse) ? mat->diffuse->value(u, v) : float3(0.f);
		specular = (mat->specular) ? mat->specular->value(u, v) : float3(0.f);
		transmission = (mat->transmission) ? mat->transmission->value(u, v) : float3(0.f);
		emission = mat->emission;
	}
};

class Intersectable {
public:
	virtual bool Intersect(const Ray& ray, Hit& hit) const = 0;

	// Surface area of the Shape
	// Mostly used for Direct light
	virtual float Area() const { return 0.f; }

	// Sample a point |S| on the shape given a reference point |ref| and
	// returns the surface normal |N| at that point
	virtual void Sample(const float3& ref, float3* S, float3* N) const {}
};

// Simple XZ plane with normal = Y
// and origin O that is used to define the plane width/height in texture coordinates
class Plane : public Intersectable {
public:
	Plane(float3 o, float2 size, shared_ptr<Material> m) : O(o), HalfSize(size / 2), mat(m) {}

	virtual bool Intersect(const Ray& ray, Hit& hit) const override {
		if (ray.D.y == 0) return false; // ray is parallel to the plane
		float t = (O.y - ray.O.y) / ray.D.y;
		if (t <= 0 || t >= ray.t) return false; // make sure intersection is close enough
		// compute intersection position
		float3 P = ray.at(t);
		// make sure intersection is inside plane
		float u = (P.x - O.x) / HalfSize.x;
		float v = (P.z - O.z) / HalfSize.y;
		if (fabsf(u) <= 1 && fabs(v) <= 1) {
			hit.u = (u + 1) / 2.0f;
			hit.v = (v + 1) / 2.0f;
			hit.mat = mat.get();
			ray.t = t;
			hit.I = P;
			hit.N = hit.ShadingN = make_float3(0, -1, 0);
			return true;
		}
		return false;
	}

	virtual float Area() const override { return 4.f * HalfSize.x * HalfSize.y; }

	virtual void Sample(const float3& ref, float3* S, float3* N) const override {
		float3 s = O;
		s.x += RandomFloat(-1, 1) * HalfSize.x;
		s.z += RandomFloat(-1, 1) * HalfSize.y;
		*S = s;
		// ensure light normal faces reference point
		float3 L = normalize(s - ref);
		float3 lightNormal = float3(0, -1, 0);
		float cos_o = dot(-L, lightNormal);
		if (cos_o < 0) lightNormal = -lightNormal;
		*N = lightNormal;
	}

public:
	float3 O;
	float2 HalfSize;
	shared_ptr<Material> mat;
};

class Sphere : public Intersectable {
public:
	Sphere(float3 center, float radius, shared_ptr<Material> m) : Center(center), r(radius), r2(radius* radius), mat(m) {}

	virtual bool Intersect(const Ray& ray, Hit& hit) const override {
		float3 oc = ray.O - Center;
		// ray.D is normalized => a = 1
		float half_b = dot(oc, ray.D);
		float c = sqrLength(oc) - r2;

		float discriminant = half_b * half_b - c;
		if (discriminant < 0) return false;
		float sqrtd = sqrtf(discriminant);

		// Find the nearest root that lies in the acceptable range.
		float root = -half_b - sqrtd;
		if (root < 0 || ray.t < root) {
			root = -half_b + sqrtd;
			if (root < 0 || ray.t < root)
				return false;
		}

		ray.t = root;
		hit.mat = mat.get();
		hit.I = ray.at(root);
		hit.ShadingN = hit.N = normalize(hit.I - Center);

		// compute point partial derivatives
		float3 pHit = hit.I - Center; // intersection point in Sphere's local coordinates system
		if (pHit.x == 0 && pHit.y == 0) pHit.x = EPSILON * r;
		float phi = atan2(pHit.y, pHit.x);
		if (phi < 0) phi += TWOPI;
		float theta = acos(clamp(pHit.z / r, -1.f, 1.f));
		float zRadius = sqrt(pHit.x * pHit.x + pHit.y * pHit.y);
		float invZRadius = 1 / zRadius;
		float cosPhi = pHit.x * invZRadius;
		float sinPhi = pHit.y * invZRadius;
		hit.dpdu = float3(-TWOPI * pHit.y, TWOPI * pHit.x, 0);
		hit.dpdv = PI * float3(pHit.z * cosPhi, pHit.z * sinPhi, -r * sin(theta));

		return true;
	}

	virtual float Area() const override { return 4.f * PI * r2; }

	virtual void Sample(const float3& ref, float3* S, float3* N) const override {
		*S = Center + RandomInSphere(r);
		*N = normalize(ref - Center);
	}

public:
	float3 Center;
	float r, r2;
	shared_ptr<Material> mat;
};
