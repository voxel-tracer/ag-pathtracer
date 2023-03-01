#pragma once

#include "microfacet.h"

class BSDF {
public:
    BSDF() {}
    BSDF(const float3& N, const float3& dpdu, MicrofacetReflection* bxdf, float eta = 1.f)
        : eta(eta),
        ng(N),
        ns(N),
        ss(normalize(dpdu)),
        ts(cross(ns, ss)),
        bxdf(bxdf) {}
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
        if (reflect) f += bxdf->f(wo, wi);
        return f;
    }
    float3 Sample_f(const float3& woWorld, float3* wiWorld, const float2& u, float* pdf) const {
        float3 wi, wo = WorldToLocal(woWorld);
        if (wo.z == 0) return 0.;
        *pdf = 0;
        float3 f = bxdf->Sample_f(wo, &wi, u, pdf);
        if (*pdf == 0) return 0;
        *wiWorld = LocalToWorld(wi);
        bool reflect = dot(*wiWorld, ng) * dot(woWorld, ng) > 0;
        if (!reflect) return float3(0.f);
        return f;
    }
    float Pdf(const float3& woWorld, const float3& wiWorld) const {
        float3 wo = WorldToLocal(woWorld), wi = WorldToLocal(wiWorld);
        if (wo.z == 0) return 0.;
        return bxdf->Pdf(wo, wi);
    }

private:
    float eta;
    float3 ng, ns, ss, ts; // shading coordinate system
    MicrofacetReflection* bxdf;
};