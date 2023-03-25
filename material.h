#pragma once

#include "disney.h"


class Material {
public:
	virtual void SetupBSDF(BSDF* bsdf) const = 0;
};

// Only support diffuse reflection for now
class DisneyMaterial : public Material {
public:
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

	void SetupBSDF(BSDF* bsdf) const {
		if (diffuse) 
			bsdf->AddBxDF(diffuse.get());
		if (retro) 
			bsdf->AddBxDF(retro.get());
		if (microfacet) 
			bsdf->AddBxDF(microfacet.get());
	}

	static std::shared_ptr<DisneyMaterial> Make(const float3& color, float roughness, float metallic) {
		return std::make_shared<DisneyMaterial>(color, roughness, metallic);
	}

private:
	float eta = 1.5f;

	shared_ptr<DisneyDiffuse> diffuse;
	shared_ptr<DisneyRetro> retro;
	shared_ptr<MicrofacetReflection> microfacet;
};

class MirrorMaterial : public Material {
public:
	MirrorMaterial(const float3& r) {
		Fresnel* fresnel = new FresnelNoOp();
		reflection = std::make_shared<SpecularReflection>(r, fresnel);
	}

	void SetupBSDF(BSDF* bsdf) const {
		bsdf->AddBxDF(reflection.get());
	}

	static std::shared_ptr<MirrorMaterial> Make(const float3& r) {
		return std::make_shared<MirrorMaterial>(r);
	}
private:
	std::shared_ptr<SpecularReflection> reflection;
};