#pragma once

class Hit {
public:
	float3 I; // intersection point
	float3 N; // normal at intersection

	const Material* mat;
	float u;
	float v;
};

class Intersectable {
public:
	virtual bool Intersect(const Ray& ray, Hit& hit) const = 0;
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
			hit.N = make_float3(0, -1, 0);
			return true;
		}
		return false;
	}

public:
	float3 O;
	float2 HalfSize;
	shared_ptr<Material> mat;
};

class Sphere : public Intersectable {
public:
	Sphere(float3 center, float radius, shared_ptr<Material> m) : Center(center), r2(radius* radius), mat(m) {}

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
		hit.N = normalize(hit.I - Center);
		return true;
	}

public:
	float3 Center;
	float r2;
	shared_ptr<Material> mat;
};
