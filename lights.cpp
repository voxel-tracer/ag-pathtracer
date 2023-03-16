#include "precomp.h"
#include "camera.h"
#include "intersectable.h"

#include "lights.h"
#include "integrator.h"

bool VisibilityTester::Unoccluded(const Scene& scene) const {
    SurfaceInteraction hit;
    return !scene.NearestIntersection(ray, hit);
}

float3 InfiniteAreaLight::Sample_Li(const SurfaceInteraction& ref, const float2& u, float3* wi, float* pdf, VisibilityTester* vis) const {
    // Sample a random direction in the hemisphere
    *wi = RandomInHemisphere(ref.shading.n);
    // Compute PDF for sampled infinite light direction
    *pdf = 1 / (2 * PI);

    // Return radiance value for infinite light direction
    *vis = VisibilityTester(Ray(ref.p + EPSILON * *wi, *wi));
    return L;
}

float InfiniteAreaLight::Pdf_Li(const SurfaceInteraction& ref, const float3& wi) const {
    return dot(ref.n, wi) > 0 ? 1 / (2 * PI) : 0.f;
}

float3 AreaLight::Sample_Li(const SurfaceInteraction& ref, const float2& u, float3* wi, float* pdf, VisibilityTester* vis) const {
    Interaction pShape = Shape->Sample(ref, u, pdf);
    if (*pdf == 0 || sqrLength(pShape.p - ref.p) == 0) {
        *pdf = 0;
        return float3(0.f);
    }
    *wi = pShape.p - ref.p;
    float dist = length(*wi);
    *wi /= dist;
    *vis = VisibilityTester(Ray(ref.p + EPSILON * *wi, *wi, dist - 10 * EPSILON)); // TODO find why we need 10x to avoid intersecting the light
    return L;
}

float AreaLight::Pdf_Li(const SurfaceInteraction& ref, const float3& wi) const {
    return Shape->Pdf(ref, wi);
}
