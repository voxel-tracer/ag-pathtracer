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

enum MaterialType { DIFFUSE, GLASS, MIRROR };

class Material {
public:
	Material() : baseColor(make_shared<SolidColor>(float3(1))), type(DIFFUSE), ref_idx(1.0f) {}
	Material(shared_ptr<Texture> c, MaterialType type = DIFFUSE, float ref_idx = 1.0f) : baseColor(c), type(type), ref_idx(ref_idx) {}

public:
	float ref_idx;	// reflection index
	shared_ptr<Texture> baseColor;	// albedo for DIFFUSE, reflection for MIRROR and absorption for GLASS
	MaterialType type;

	static shared_ptr<Material> make_glass(float ref_idx) {
		return make_shared<Material>(make_shared<SolidColor>(0.f), GLASS, ref_idx);
	}

	static shared_ptr<Material> make_mirror(float r, float g, float b) {
		return make_shared<Material>(make_shared<SolidColor>(r, g, b), MIRROR);
	}
};
