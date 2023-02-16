#pragma once

class Light {
public:
    Light(float3 l) : L(l) {}

    virtual bool Sample(const float3& ref, float3* S, float3* N) const { return false; }

    virtual float Area() const { return 0.f; }

    virtual float3 Le(const Ray& ray) const {
        return float3(0.f);
    }

    float3 L;
};

class InfiniteAreaLight : public Light {
public:
    InfiniteAreaLight(float3 l) : Light(l) {}

    virtual float3 Le(const Ray& ray) const override {
        return L;
    }
};

class PointLight : public Light {
public:
    PointLight(const float3& pos, const float3& l) : Light(l), Pos(pos) {}

    virtual bool Sample(const float3& ref, float3* S, float3* N) const override {
        *S = Pos;
        float3 L = normalize(Pos - ref);
        *N = -L;
        return true;
    }

    float3 Pos;
};

class AreaLight : public Light {
public:
    AreaLight(shared_ptr<Intersectable> shape, const float3& l) : Light(l), Shape(shape) {}

    virtual bool Sample(const float3& ref, float3* S, float3* N) const override {
        Shape->Sample(ref, S, N);
        return true;
    }

    virtual float Area() const override { return Shape->Area(); }

    shared_ptr<Intersectable> Shape;
};