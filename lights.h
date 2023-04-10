#pragma once

class Scene;
class Ray;
class SurfaceInteraction;
class HDRTexture;
struct Distribution1D;

#define ILS

class VisibilityTester {
public:
    VisibilityTester(): ray() {}
    VisibilityTester(const Ray& ray) :ray(ray) {}

    bool Unoccluded(const Scene& scene) const;

private:
    Ray ray;
};

class Light {
public:
    virtual float3 Sample_Li(const SurfaceInteraction& ref, const float2& u, float3* wi, float* pdf, VisibilityTester* vis) const {
        return float3(0.f);
    }

    virtual float Pdf_Li(const SurfaceInteraction& si, const float3& wi) const {
        return 0.f;
    }

    virtual float3 Le(const Ray& ray) const { return float3(0.f); }

    virtual bool IsInfinite() const { return false; }
};

class UniformInfiniteLight : public Light {
public:
    UniformInfiniteLight(const float3& l) : Lemit(l) {}

    float3 Sample_Li(const SurfaceInteraction& ref, const float2& u, float3* wi, float* pdf, VisibilityTester* vis) const;

    float Pdf_Li(const SurfaceInteraction& si, const float3& wi) const;

    float3 Le(const Ray& ray) const { return Lemit; }

    virtual bool IsInfinite() const { return true; }

protected:
    const float3 Lemit;
};

class InfiniteAreaLight : public Light {
public:
    InfiniteAreaLight(const std::string& texmap);

    float3 Sample_Li(const SurfaceInteraction& ref, const float2& u, float3* wi, float* pdf, VisibilityTester* vis) const;

    float Pdf_Li(const SurfaceInteraction& si, const float3& wi) const;

    float3 Le(const Ray& ray) const;

    virtual bool IsInfinite() const { return true; }

private:
    std::shared_ptr<HDRTexture> Lmap;
#ifdef ILS
    std::shared_ptr<Distribution1D> distrib;
#endif
};

class AreaLight : public Light {
public:
    AreaLight(shared_ptr<Intersectable> shape, const float3& l) : Lemit(l), Shape(shape) {
        shape->SetAreaLight(this);
    }

    float3 Sample_Li(const SurfaceInteraction& ref, const float2& u, float3* wi, float* pdf, VisibilityTester* vis) const;

    float Pdf_Li(const SurfaceInteraction& si, const float3& wi) const;

    float3 L(const SurfaceInteraction& si, const float3& w) const { return Lemit; }

    shared_ptr<Intersectable> Shape;
protected:
    const float3 Lemit;
};