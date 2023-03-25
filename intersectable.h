#pragma once

#include "material.h"
#include "camera.h"
#include "reflection.h"

class BSDF;

class SurfaceInteraction;
class AreaLight;

struct Interaction {
	float3 p;
	float3 n;
};

class Intersectable {
public:
	Intersectable(shared_ptr<Material> material) : material(material) {}

	virtual bool Intersect(const Ray& ray, SurfaceInteraction& hit) const = 0;

	// Surface area of the Shape
	// Mostly used for Direct light
	virtual float Area() const { return 0.f; }

	// Sample a point uniformly on the surface of the shape
	virtual Interaction Sample(const float2& u, float* pdf) const {
		*pdf = 0;
		return Interaction();
	}

	// Sample a point on the shape surface according to a probability density
	// with respect to solid angle form the reference point _ref_
	virtual Interaction Sample(const SurfaceInteraction& ref, const float2& u, float* pdf) const {
		*pdf = 0;
		return Interaction();
	}

	virtual float Pdf(const SurfaceInteraction& ref, const float3& wi) const {
		return 0.f;
	}

	AreaLight* GetAreaLight() const {
		return arealight.get();
	}

	void SetAreaLight(shared_ptr<AreaLight> light) {
		arealight = light;
	}

	std::shared_ptr<Material> material;
	std::shared_ptr<AreaLight> arealight;
};

class SurfaceInteraction {
public:
	SurfaceInteraction() {}
	SurfaceInteraction(const float3& p, const float2& uv, const float3& wo, const float3& dpdu, const float3& dpdv, const Intersectable* shape) : 
			p(p), 
			wo(wo), 
			n(normalize(cross(dpdu, dpdv))),
			uv(uv), 
			dpdu(dpdu),
			dpdv(dpdv),
			shape(shape) {
		// initialize the shading geometry
		shading.n = n;
		shading.dpdu = dpdu;
		shading.dpdv = dpdv;
	}

	void SetShadingGeometry(const float3& dpdus, const float3& dpdvs, bool orientationIsAuthorative) {
		// Compute shading normal for surface interaction
		shading.n = normalize(cross(dpdus, dpdvs));
		if (orientationIsAuthorative)
			n = Faceforward(n, shading.n);
		else
			shading.n = Faceforward(shading.n, n);
		dpdu = dpdus;
		dpdv = dpdvs;
	}

	void EvalMaterial() {
		Material* mat = shape->material.get();

		if (mat->disney) {
			bsdf = BSDF(*this);
			if (mat->disney->diffuse)
				bsdf.AddBxDF(mat->disney->diffuse.get());
			if (mat->disney->retro)
				bsdf.AddBxDF(mat->disney->retro.get());
			if (mat->disney->microfacet)
				bsdf.AddBxDF(mat->disney->microfacet.get());
			hasBSDF = true;
		}
		else {
			hasBSDF = false;
		}

		emission = mat->emission;
	}

	// Interaction public data
	float3 p; // intersection point
	float3 wo;
	float3 n;

	float2 uv;
	float3 dpdu, dpdv;
	const Intersectable* shape = nullptr;

	struct {
		float3 n;
		float3 dpdu, dpdv;
	} shading;

	BSDF bsdf;
	bool hasBSDF;
	float3 emission;
};

// Simple XZ plane with normal = Y
// and origin O that is used to define the plane width/height in texture coordinates
class Plane : public Intersectable {
public:
	Plane(float3 o, float2 size, shared_ptr<Material> m) : Intersectable(m), O(o), HalfSize(size / 2) {}

	virtual bool Intersect(const Ray& ray, SurfaceInteraction& hit) const override {
		if (ray.D.y == 0) return false; // ray is parallel to the plane
		float t = (O.y - ray.O.y) / ray.D.y;
		if (t <= 0 || t >= ray.t) return false; // make sure intersection is close enough
		// compute intersection position
		float3 P = ray.at(t);
		// make sure intersection is inside plane
		float u = (P.x - O.x) / HalfSize.x;
		float v = (P.z - O.z) / HalfSize.y;
		if (fabsf(u) <= 1 && fabs(v) <= 1) {
			hit = SurfaceInteraction(P, float2((u + 1) * .5f, (v + 1) * .5f), -ray.D, float3(0, 0, 1), float3(1, 0, 0), this);
			ray.t = t;
			return true;
		}
		return false;
	}

	virtual float Area() const override { return 4.f * HalfSize.x * HalfSize.y; }

public:
	float3 O;
	float2 HalfSize;
};

class Sphere : public Intersectable {
public:
	Sphere(float3 center, float radius, shared_ptr<Material> m) : 
		Intersectable(m), Center(center), r(radius), r2(radius* radius) {}

	virtual bool Intersect(const Ray& ray, SurfaceInteraction& hit) const override {
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

		// compute point partial derivatives
		float3 p = ray.at(root);
		float3 pHit = p - Center; // intersection point in Sphere's local coordinates system
		if (pHit.x == 0 && pHit.y == 0) pHit.x = EPSILON * r;
		float phi = std::atan2(pHit.y, pHit.x);
		if (phi < 0) phi += TWOPI;
		float u = phi * INV2PI;

		float theta = std::acos(clamp(pHit.z / r, -1.f, 1.f));
		float v = theta * INVPI;

		float zRadius = sqrt(pHit.x * pHit.x + pHit.y * pHit.y);
		float invZRadius = 1 / zRadius;
		float cosPhi = pHit.x * invZRadius;
		float sinPhi = pHit.y * invZRadius;
		float3 dpdu = float3(-TWOPI * pHit.y, TWOPI * pHit.x, 0);
		float3 dpdv = PI * float3(pHit.z * cosPhi, pHit.z * sinPhi, -r * sin(theta));

		// NOTE: I don't know why but I need to invert the computed dpdu/dpdv in order to get the normal to correctly point out of the sphere
		hit = SurfaceInteraction(p, float2(u, v), -ray.D, dpdv, dpdu, this);
		ray.t = root;

		return true;
	}

	virtual float Area() const override { return 4.f * PI * r2; }

	virtual Interaction Sample(const float2& u, float* pdf) const override {
		float3 pObj = Center + r * RandomInSphere(u);
		Interaction it;
		it.n = normalize(pObj);
		it.p = pObj;
		*pdf = 1 / Area();
		return it;
	}

	virtual Interaction Sample(const SurfaceInteraction& ref, const float2& u, float* pdf) const override {
		float3 pCenter = Center;

		// Sample uniformly on sphere if pt is inside it
		float3 pOrigin = ref.p;
		if (sqrLength(pOrigin - pCenter) <= r2) {
			Interaction intr = Sample(u, pdf);
			float3 wi = intr.p - ref.p;
			if (sqrLength(wi) == 0)
				*pdf = 0;
			else {
				// convert from area measure returned by Sample() call above to
				// solid angle measure
				wi = normalize(wi);
				*pdf *= sqrLength(ref.p - intr.p) / absdot(intr.n, -wi);
			}
			if (std::isinf(*pdf))*pdf = 0;
			return intr;
		}

		// Sample sphere uniformly inside subtended cone

		// Compute coordinate system for sphere sampling
		float dc = length(ref.p - pCenter);
		float invDc = 1 / dc;
		float3 wc = (pCenter - ref.p) * invDc;
		float3 wcX, wcY;
		CoordinateSystem(wc, &wcX, &wcY);

		// Compute _theta_ and _phi_ values for sample in cone
		float sinThetaMax = r * invDc;
		float sinThetaMax2 = sinThetaMax * sinThetaMax;
		float invSinThetaMax = 1 / sinThetaMax;
		float cosThetaMax = std::sqrt(std::max(0.f, 1 - sinThetaMax2));

		float cosTheta = (cosThetaMax - 1) * u[0] + 1;
		float sinTheta2 = 1 - cosTheta * cosTheta;

		if (sinThetaMax2 < 0.00068523f /* sin^2(1.5 deg) */) {
			/* Fall back to a Taylor series expansion for small angles, where
			   the standard approach suffers from severe cancellation errors */
			sinTheta2 = sinThetaMax2 * u[0];
			cosTheta = std::sqrt(1 - sinTheta2);
		}

		// Compute angle _alpha_ from center of sphere to sampled point on surface
		float cosAlpha = sinTheta2 * invSinThetaMax +
			cosTheta * std::sqrt(std::max(0.f, 1.f - sinTheta2 * invSinThetaMax * invSinThetaMax));
		float sinAlpha = std::sqrt(std::max(0.f, 1.f - cosAlpha * cosAlpha));
		float phi = u[1] * 2 * PI;

		// Compute surface normal and sampled point on sphere
		float3 nWorld =
			SphericalDirection(sinAlpha, cosAlpha, phi, -wcX, -wcY, -wc);
		float3 pWorld = pCenter + r * float3(nWorld.x, nWorld.y, nWorld.z);

		// Return _Interaction_ for sampled point on sphere
		Interaction it;
		it.p = pWorld;
		it.n = nWorld;

		// Uniform cone PDF.
		*pdf = 1 / (2 * PI * (1 - cosThetaMax));

		return it;
	}

	virtual float Pdf(const SurfaceInteraction& ref, const float3& wi) const override {
		float3 pCenter = Center;
		// Return uniform PDF if point is inside the Sphere
		float3 pOrigin = ref.p;
		if (sqrLength(pOrigin - pCenter) <= r2)
			return 1 / (4 * PI); // oversimplification as we don't really care about this case for now

		// Compute general sphere PDF
		float sinThetaMax2 = r2 / sqrLength(ref.p - pCenter);
		float costThetaMax = std::sqrt(std::max(0.f, 1 - sinThetaMax2));
		return UniformConePdf(costThetaMax);
	}

public:
	float3 Center;
	float r, r2;
};
