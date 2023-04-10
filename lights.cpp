#include "precomp.h"
#include "camera.h"
#include "intersectable.h"

#include "lights.h"
#include "integrator.h"
#include "texture.h"
#include "sampling.h"

bool VisibilityTester::Unoccluded(const Scene& scene) const {
    return !scene.IntersectP(ray);
}


float3 UniformInfiniteLight::Sample_Li(const SurfaceInteraction& ref, const float2& u, float3* wi, float* pdf, VisibilityTester* vis) const {
    // Sample a random direction in the hemisphere
    *wi = RandomInHemisphere(ref.shading.n);
    // Compute PDF for sampled infinite light direction
    *pdf = INV2PI;

    // Return radiance value for infinite light direction
    *vis = VisibilityTester(Ray(ref.p + EPSILON * *wi, *wi));
    return Lemit;
}

float UniformInfiniteLight::Pdf_Li(const SurfaceInteraction& ref, const float3& wi) const {
    return dot(ref.n, wi) > 0 ? INV2PI : 0.f;
}


InfiniteAreaLight::InfiniteAreaLight(const std::string& texmap) {
    Lmap = std::make_shared<HDRTexture>(texmap);
#ifdef ILS
    // Compute light map scaled by sin(theta)
    std::vector<float> pdf(Lmap->Width() * Lmap->Height());
    for (int idx = 0; idx < pdf.size(); idx++) {
        int x = idx % Lmap->Width();
        int y = idx / Lmap->Width();
        float th = (y + .5f) * PI / Lmap->Height();
        float3 value = Lmap->GetPixel(x, y);
        float maxComponent = max(value.x, max(value.y, value.z));
        // multiply pixel value by sin(th) to account for the row stretch near the pole
        pdf[idx] = maxComponent * std::sin(th);
    }
    // Initialize 1D distribution
    distrib = std::make_shared<Distribution1D>(&pdf[0], (int)pdf.size());
#endif
}

float3 InfiniteAreaLight::Sample_Li(const SurfaceInteraction& ref, const float2& u, float3* wi, float* pdf, VisibilityTester* vis) const {
#ifdef ILS
    // Use Distribution1D to sample a direction in the sphere
    float mapPdf;
    float sample = distrib->SampleContinuous(RandomFloat(), &mapPdf);
    if (mapPdf == 0) return float3(0.f);

    // compute index to texture coordinates
    int idx = int(sample * distrib->Count());
    float2 uv(
        ((idx % Lmap->Width()) + .5f) / Lmap->Width(),
        ((idx / Lmap->Width()) + .5f) / Lmap->Height()
    );
    // Convert infinite light sample point to direction
    float theta = uv.y * PI, phi = uv.x * TWOPI;
    float cosTheta = std::cos(theta), sinTheta = std::sin(theta);
    float sinPhi = std::sin(phi), cosPhi = std::cos(phi);
    *wi = float3(sinTheta*cosPhi, cosTheta, sinTheta*sinPhi);

    // Compute PDF for sampled infinite light direction
    *pdf = mapPdf / (2 * PI * PI * sinTheta);
    if (sinTheta == 0) *pdf = 0;

    // Return radiance value for infinite light direction
    Ray ray(ref.p + EPSILON * ref.n, *wi);
    *vis = VisibilityTester(ray);

    return Le(ray);
#else
    // Sample a random direction in the hemisphere
    * wi = RandomInHemisphere(ref.shading.n);
    // Compute PDF for sampled infinite light direction
    *pdf = INV2PI;

    // Return radiance value for infinite light direction
    Ray ray(ref.p + EPSILON * ref.n, *wi);
    *vis = VisibilityTester(ray);
    return Le(ray);

#endif
}

float InfiniteAreaLight::Pdf_Li(const SurfaceInteraction& si, const float3& wi) const {
#ifdef ILS
    // Find pixel coordinates corresponding to direction _wi_
    float3 w = normalize(wi); w = float3(w.x, w.z, w.y);
    float theta = SphericalTheta(w), phi = SphericalPhi(w);
    float sinTheta = std::sin(theta);
    if (sinTheta == 0) return 0;
    int x = clamp(int(phi * INV2PI * Lmap->Width()), 0, Lmap->Width() - 1);
    int y = clamp(int(theta * INVPI * Lmap->Height()), 0, Lmap->Height() - 1);

    return distrib->Count() * distrib->DiscretePDF(y * Lmap->Width() + x) / (2 * PI * PI * sinTheta);
#else
    return dot(si.n, wi) > 0 ? INV2PI : 0.f;
#endif
}

float3 InfiniteAreaLight::Le(const Ray& ray) const {
    float3 w = normalize(ray.D); w = float3(w.x, w.z, w.y);
    float2 st(SphericalPhi(w) * INV2PI, SphericalTheta(w) * INVPI);
    return Lmap->value(st.x, st.y);
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
    return Lemit;
}

float AreaLight::Pdf_Li(const SurfaceInteraction& ref, const float3& wi) const {
    return Shape->Pdf(ref, wi);
}
