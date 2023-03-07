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

// Only support diffuse reflection for now
struct DisneyMaterial {
	float eta = 1.5f;

	shared_ptr<DisneyDiffuse> diffuse;
	shared_ptr<DisneyRetro> retro;
	shared_ptr<MicrofacetReflection> microfacet;

	DisneyMaterial(const float3& color, float roughness, float metallic) {
		// Diffuse
		float3 c = color;
		float metallicWeight = metallic;
		float e = eta;
		float strans = 0.f;
		float diffuseWeight = (1 - metallicWeight) * (1 - strans);
		float dt = 0.f; // 0: all diffuse is reflected -> 1, transmitted
		float rough = roughness; 

		//TODO we need to properly compute the following
		float3 Ctint = float3(1.f);

		if (diffuseWeight > 0) {
			float3 sd(0.f);
			if (IsBlack(sd)) {
				// No subsurface scattering; use regular (Fresnel modified) diffuse
				diffuse = make_shared<DisneyDiffuse>(diffuseWeight * c);
			}

			// Retro-reflection
			retro = make_shared<DisneyRetro>(diffuseWeight * c, rough);
		}

		// Create the microfacet distribution for metallic and/or specular transmission
		float aspect = 1.f;
		float ax = std::max(.001f, sqr(rough) / aspect);
		float ay = std::max(.001f, sqr(rough) * aspect);
		MicrofacetDistribution* distrib = new DisneyMicrofacetDistribution(ax, ay);

		// Specular is Trowbridge-Reitz with a modified Fresnel function.
		float specTint = 0.f;
		float3 Cspec0 = Lerp(metallicWeight, SchlickR0FromEta(e) * Lerp(specTint, float3(1.f), Ctint), c);
		Fresnel* fresnel = new DisneyFresnel(Cspec0, metallicWeight, e);
		microfacet = make_shared<MicrofacetReflection>(float3(1.), distrib, fresnel);
	}
};

class Material {
public:
	float ref_idx = 1.f;				// reflection index
	shared_ptr<Texture> diffuse;
	shared_ptr<Texture> specular;
	shared_ptr<Texture> transmission;
	shared_ptr<MicrofacetReflection> microfacet;
	shared_ptr<DisneyMaterial> disney;
	float3 emission = float3(0.f);

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

	static shared_ptr<Material> make_metal(float roughness, const float3& eta, const float3& k) {
		auto mat = make_shared<Material>();

		float uRough = TrowbridgeReitzDistribution::RoughnessToAlpha(roughness);
		float vRough = TrowbridgeReitzDistribution::RoughnessToAlpha(roughness);
		MicrofacetDistribution* distribution = new TrowbridgeReitzDistribution(uRough, vRough);
		Fresnel* fresnel = new FresnelConductor(float3(1.f), eta, k);
		mat->microfacet = make_shared<MicrofacetReflection>(float3(1.f), distribution, fresnel);
		return mat;
	}

	static shared_ptr<Material> make_disney(const float3& color, float roughness, float metallic) {
		auto mat = make_shared<Material>();
		mat->disney = make_shared<DisneyMaterial>(color, roughness, metallic);
		return mat;
	}
};
