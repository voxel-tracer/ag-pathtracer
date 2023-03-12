#include "precomp.h"

#include "reflection.h"
#include "intersectable.h"

BSDF::BSDF(const SurfaceInteraction& si, float eta) : 
    eta(eta),
    ng(si.n),
    ns(si.shading.n),
    ss(normalize(si.shading.dpdu)),
    ts(cross(ns, ss)) {}
