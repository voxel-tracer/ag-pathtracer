#pragma once

class Texture {
public:
	virtual float3 value(float u, float v) const = 0;
};

class SolidColor : public Texture {
public:
	SolidColor(float a) : SolidColor(a, a, a) {}
	SolidColor(float r, float g, float b) : SolidColor(float3(r, g, b)) {}
	SolidColor(float3 c) : color(c) {}

	virtual float3 value(float u, float v) const override {
		return color;
	}

	static shared_ptr<SolidColor> make(float3 v) {
		return make_shared<SolidColor>(v);
	}

private:
	float3 color;
};

class CheckerTexture : public Texture {
public:
	CheckerTexture(shared_ptr<Texture> c1, shared_ptr<Texture> c2) : color1(c1), color2(c2) {}

	virtual float3 value(float u, float v) const override {
		if (((int)(u * 8) + (int)(v * 8)) % 2)
			return color1->value(u, v);
		return color2->value(u, v);
	}

private:
	shared_ptr<Texture> color1;
	shared_ptr<Texture> color2;
};

class Material {
public:
	float ref_idx = 1.f;				// reflection index
	shared_ptr<Texture> diffuse;
	shared_ptr<Texture> specular;
	shared_ptr<Texture> transmission;
	float3 emission;

	static shared_ptr<Material> make_lambertian(shared_ptr<Texture> diffuse) {
		auto lambertian = make_shared<Material>();
		lambertian->diffuse = diffuse;
		return lambertian;
	}

	static shared_ptr<Material> make_glass(float ref_idx, float3 transmission = float3(1.f)) {
		auto glass = make_shared<Material>();
		glass->ref_idx = ref_idx;
		glass->transmission = SolidColor::make(transmission);
		glass->specular = SolidColor::make(1.0f);
		return glass;
	}

	static shared_ptr<Material> make_mirror(float r, float g, float b) {
		auto mirror = make_shared<Material>();
		mirror->specular = SolidColor::make(float3(r, g, b));
		return mirror;
	}

	static shared_ptr<Material> make_emitter(float r, float g, float b) {
		auto mat = make_shared <Material>();
		mat->emission = float3(r, g, b);
		return mat;
	}
};
