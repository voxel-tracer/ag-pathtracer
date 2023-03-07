#pragma once

#include "reflection.h"

// https://seblagarde.wordpress.com/2013/04/29/memo-on-fresnel-equations/
//
// The Schlick Fresnel approximation is:
//
// R = R(0) + (1 - R(0)) (1 - cos theta)^5,
//
// where R(0) is the reflectance at normal indicence.
inline float SchlickWeight(float cosTheta) {
    float m = clamp(1 - cosTheta, 0.f, 1.f);
    return (m * m) * (m * m) * m;
}

inline float3 FrSchlick(const float3& R0, float cosTheta) {
    return Lerp(SchlickWeight(cosTheta), R0, float3(1.f));
}

// For a dielectric, R(0) = (eta - 1)^2 / (eta + 1)^2, assuming we're
// coming from air..
inline float SchlickR0FromEta(float eta) { return sqr(eta - 1) / sqr(eta + 1); }

class DisneyDiffuse : public BxDF {
public:
    DisneyDiffuse(const float3& R) : R(R) {}
    virtual float3 f(const float3& wo, const float3& wi) const override {
        float Fo = SchlickWeight(AbsCosTheta(wo)),
            Fi = SchlickWeight(AbsCosTheta(wi));

        // Diffuse fresnel - go from 1 at normal incidence to .5 at grazing.
        // Burley 2015, eq (4).
        return R * INVPI * (1 - Fo / 2) * (1 - Fi / 2);
    }
private:
    float3 R;
};

class DisneyRetro : public BxDF {
public:
    DisneyRetro(const float3& R, float roughness) : R(R), roughness(roughness) {}
    virtual float3 f(const float3& wo, const float3& wi) const override {
        float3 wh = wi + wo;
        if (wh.x == 0 && wh.y == 0 && wh.z == 0) return float3(0.f);
        wh = normalize(wh);
        float cosThetaD = dot(wi, wh);

        float Fo = SchlickWeight(AbsCosTheta(wo)),
            Fi = SchlickWeight(AbsCosTheta(wi));
        float Rr = 2 * roughness * cosThetaD * cosThetaD;

        // Burley 2015, eq (4).
        return R * INVPI * Rr * (Fo + Fi + Fo * Fi * (Rr - 1));
    }

private:
    float3 R;
    float roughness;
};

class DisneyFresnel : public Fresnel {
public:
    DisneyFresnel(const float3& R0, float metallic, float eta) : R0(R0), metallic(metallic), eta(eta) {}
    virtual float3 Evaluate(float cosI) const override {
        return Lerp(metallic, float3(FrDielectric(cosI, 1, eta)), FrSchlick(R0, cosI));
    }
private:
    const float3 R0;
    const float metallic, eta;
};

class DisneyMicrofacetDistribution : public TrowbridgeReitzDistribution {
public:
    DisneyMicrofacetDistribution(float alphax, float alphay)
        : TrowbridgeReitzDistribution(alphax, alphay) {}

    virtual float G(const float3& wo, const float3& wi) const override {
        // Disney uses the separable masking-shadowing model.
        return G1(wo) * G1(wi);
    }
};