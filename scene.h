#pragma once

class Scene {
public:
	bool Intersect(const Ray& ray, SurfaceInteraction& hit) const {
		bool found = false;
		for (auto& primitive : primitives) {
			if (primitive->Intersect(ray, hit)) {
				found = true;
			}
		}
		return found;
	}

	bool IntersectP(const Ray& ray) const {
		for (auto& primitive : primitives)
			if (primitive->IntersectP(ray)) return true;
		return false;
	}

	void addAreaLight(std::shared_ptr<Intersectable> shape, const float3& L) {
		primitives.push_back(shape);
		auto area = make_shared<AreaLight>(shape, L);
		lights.push_back(area);
	}

	vector<shared_ptr<Intersectable>> primitives;
	vector<shared_ptr<Light>> lights;
	CameraDesc camera;
};
