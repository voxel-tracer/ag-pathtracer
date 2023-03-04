#pragma once

#include "microfacet.h"

class BxDF {
public:
    virtual float3 f(const float3& wo, const float3& wi) const = 0;
    virtual float3 Sample_f(const float3& wo, float3* wi, const float2& u, float* pdf) const {
        // Cosine-sample the hemisphere, flipping the direction if necessary
        //*wi = CosineWeightedRandomInHemisphere(float3(0, 0, 1));
        *wi = CosineSampleHemisphere(u);
        if (wo.z < 0) wi->z *= -1;
        *pdf = Pdf(wo, *wi);
        return f(wo, *wi);
    }
    virtual float Pdf(const float3& wo, const float3& wi) const {
        return SameHemisphere(wo, wi) ? AbsCosTheta(wi) * INVPI : 0;
    }
};

class MicrofacetReflection : public BxDF {
public:
    MicrofacetReflection(const float3& R, MicrofacetDistribution* distribution, Fresnel* fresnel)
        : R(R), distribution(distribution), fresnel(fresnel) {}
    virtual float3 f(const float3& wo, const float3& wi) const override {
        float cosThetaO = AbsCosTheta(wo), cosThetaI = AbsCosTheta(wi);
        float3 wh = wi + wo;
        // Handle degenerate cases for microfacet reflection
        if (cosThetaI == 0 || cosThetaO == 0) return float3(0.f);
        if (wh.x == 0 && wh.y == 0 && wh.z == 0) return float3(0.f);
        wh = normalize(wh);
        // For the Fresnel call, make sure that wh is in the same hemisphere
        // as the surface normal, so that TIR is handled correctly.
        float3 F = fresnel->Evaluate(dot(wi, Faceforward(wh, float3(0, 0, 1))));
        return R * distribution->D(wh) * distribution->G(wo, wi) * F /
            (4 * cosThetaI * cosThetaO);
    }
    virtual float3 Sample_f(const float3& wo, float3* wi, const float2& u, float* pdf) const override {
        // Sample microfacet orientation $\wh$ and reflected direction $\wi$
        if (wo.z == 0) return 0.;
        float3 wh = distribution->Sample_wh(wo, u);
        if (dot(wo, wh) < 0) return 0.;   // Should be rare
        *wi = Reflect(wo, wh);
        if (!SameHemisphere(wo, *wi)) return float3(0.f);

        // Compute PDF of _wi_ for microfacet reflection
        *pdf = distribution->Pdf(wo, wh) / (4 * dot(wo, wh));
        return f(wo, *wi);
    }
    virtual float Pdf(const float3& wo, const float3& wi) const override {
        if (!SameHemisphere(wo, wi)) return 0;
        float3 wh = normalize(wo + wi);
        return distribution->Pdf(wo, wh) / (4 * dot(wo, wh));
    }

private:
    const float3 R;
    // TODO simplify distribution attribute as we only support a single type of microfacet distribution
    const MicrofacetDistribution* distribution;
    const Fresnel* fresnel;
};

// Slightly simplified version that contains a fixed number of non-specular reflection BxDFs
class BSDF {
public:
    BSDF() {}
    BSDF(const float3& N, const float3& dpdu, float eta = 1.f)
        : eta(eta),
        ng(N),
        ns(N),
        ss(normalize(dpdu)),
        ts(cross(ns, ss)) {}

    void AddBxDF(BxDF* b) {
        bxdfs[nBxDFs++] = b;
    }
    int NumComponents() const { return nBxDFs; }

    float3 WorldToLocal(const float3& v) const {
        return float3(dot(v, ss), dot(v, ts), dot(v, ns));
    }
    float3 LocalToWorld(const float3& v) const {
        return float3(
            ss.x * v.x + ts.x * v.y + ns.x * v.z,
            ss.y * v.x + ts.y * v.y + ns.y * v.z,
            ss.z * v.x + ts.z * v.y + ns.z * v.z);
    }

    float3 f(const float3& woW, const float3& wiW) const {
        float3 wi = WorldToLocal(wiW), wo = WorldToLocal(woW);
        if (wo.z == 0) return 0.;
        bool reflect = dot(wiW, ng) * dot(woW, ng) > 0;
        float3 f(0.f);
        for (int i = 0; i < nBxDFs; i++)
            if (reflect) 
                f += bxdfs[i]->f(wo, wi);
        return f;
    }
    float3 Sample_f(const float3& woWorld, float3* wiWorld, const float2& u, float* pdf) const {
        // choose which BxDF to sample
        int matchingComps = NumComponents();
        if (matchingComps == 0) {
            *pdf = 0;
            return float3(0.f);
        }
        int comp = 
            std::min((int)std::floor(u[0] * matchingComps), matchingComps - 1);

        // get BxDF pointer to choosen component
        BxDF* bxdf = nullptr;
        int count = comp;
        for (int i = 0; i < nBxDFs; i++)
            if (count-- == 0) {
                bxdf = bxdfs[i];
                break;
            }

        // Remap BxDF sample u to [0,1)
        float2 uRemapped(std::min(u[0] * matchingComps - comp, FloatOneMinusEpsilon), u[1]);

        // Sample chosen BxDF
        float3 wi, wo = WorldToLocal(woWorld);
        if (wo.z == 0) return 0.;
        *pdf = 0;
        float3 f = bxdf->Sample_f(wo, &wi, uRemapped, pdf);
        if (*pdf == 0) return 0;
        *wiWorld = LocalToWorld(wi);

        // Compute overall PDF with all matching BxDFs
        if (matchingComps > 1) {
            for (int i = 0; i < nBxDFs; i++)
                if (bxdfs[i] != bxdf)
                    *pdf += bxdfs[i]->Pdf(wo, wi);
        }
        if (matchingComps > 1) *pdf /= matchingComps;

        // Compute value of BSDF for sampled direction
        bool reflect = dot(*wiWorld, ng) * dot(woWorld, ng) > 0;
        f = 0.f;
        for (int i = 0; i < nBxDFs; i++)
            if (reflect) 
                f += bxdfs[i]->f(wo, wi);
        return f;
    }

    float Pdf(const float3& woWorld, const float3& wiWorld) const {
        if (nBxDFs == 0) return 0.f;
        float3 wo = WorldToLocal(woWorld), wi = WorldToLocal(wiWorld);
        if (wo.z == 0) return 0.;
        float pdf = 0.f;
        int matchingComps = 0;
        for (int i = 0; i < nBxDFs; i++) {
            ++matchingComps;
            pdf += bxdfs[i]->Pdf(wo, wi);
        }
        float v = matchingComps > 0 ? pdf / matchingComps : 0.f;
        return v;
    }

private:
    float eta;
    float3 ng, ns, ss, ts; // shading coordinate system

    int nBxDFs = 0;
    static constexpr int MaxBxDFs = 3;
    BxDF* bxdfs[MaxBxDFs];
};