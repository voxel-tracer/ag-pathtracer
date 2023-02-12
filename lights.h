#pragma once

class Light {
public:
    virtual float3 Sample_Li(const float3& p, float3* wi, Ray* shadowRay) const = 0;

    virtual float3 Le(const Ray& ray) const {
        return float3(0.f);
    }
};

class InfiniteAreaLight : public Light {
public:
    InfiniteAreaLight(float3 l) : L(l) {}

    virtual float3 Sample_Li(const float3& p, float3* wi, Ray* shadowRay) const override {
        return float3(0.f);
    }

    virtual float3 Le(const Ray& ray) const override {
        return L;
    }

    float3 L;
};

class PointLight : public Light {
public:
    PointLight(const float3& pos, const float3& l) : Pos(pos), L(l) {}

    virtual float3 Sample_Li(const float3& p, float3* wi, Ray* shadowRay) const override {
        float3 PtoL = Pos - p;
        float dist = length(PtoL);
        
        *wi = PtoL / dist;
        *shadowRay = Ray(p + EPSILON * (*wi), *wi, dist);

        return L / (dist * dist);
    }


    float3 Pos;
    float3 L;
};