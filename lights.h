#pragma once

class Scene;
class Ray;
class SurfaceInteraction;

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
    Light(float3 l) : L(l) {}

    virtual float3 Sample_Li(const SurfaceInteraction& ref, const float2& u, float3* wi, float* pdf, VisibilityTester* vis) const {
        return float3(0.f);
    }

    virtual float Pdf_Li(const SurfaceInteraction& si, const float3& wi) const {
        return 0.f;
    }

    virtual float3 Le(const Ray& ray) const {
        return float3(0.f);
    }

    float3 L;
};

class InfiniteAreaLight : public Light {
public:
    InfiniteAreaLight(float3 l) : Light(l) {}

    virtual float3 Sample_Li(const SurfaceInteraction& ref, const float2& u, float3* wi, float* pdf, VisibilityTester* vis) const override;

    virtual float Pdf_Li(const SurfaceInteraction& si, const float3& wi) const override;

    virtual float3 Le(const Ray& ray) const override {
        return L;
    }

};

class AreaLight : public Light {
public:
    AreaLight(shared_ptr<Intersectable> shape, const float3& l) : Light(l), Shape(shape) {}

    virtual float3 Sample_Li(const SurfaceInteraction& ref, const float2& u, float3* wi, float* pdf, VisibilityTester* vis) const override;

    virtual float Pdf_Li(const SurfaceInteraction& si, const float3& wi) const override;

    shared_ptr<Intersectable> Shape;
};