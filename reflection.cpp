#include "precomp.h"

#include "reflection.h"
#include "intersectable.h"

BSDF::BSDF(const SurfaceInteraction& si, float eta) : 
    eta(eta),
    ng(si.n),
    ns(si.shading.n),
    ss(normalize(si.shading.dpdu)),
    ts(cross(ns, ss)) {}

float3 SpecularReflection::Sample_f(const float3& wo, float3* wi, const float2& u, float* pdf) const {
    // Compute perfect specular reflection direction
    *wi = float3(-wo.x, -wo.y, wo.z);
    *pdf = 1;
    return fresnel->Evaluate(CosTheta(*wi)) * R / AbsCosTheta(*wi);
}