#pragma once

// BSDF Inline Functions
// These functions assume w is in the shading coordinate space where N is (0, 0, 1)
inline float CosTheta(const float3& w) { return w.z; }
inline float Cos2Theta(const float3& w) { return w.z * w.z; }
inline float AbsCosTheta(const float3& w) { return abs(w.z); }
inline float Sin2Theta(const float3& w) {
    return max(0.f, 1.f - Cos2Theta(w));
}

inline float SinTheta(const float3& w) { return sqrt(Sin2Theta(w)); }

inline float TanTheta(const float3& w) { return SinTheta(w) / CosTheta(w); }

inline float Tan2Theta(const float3& w) {
    return Sin2Theta(w) / Cos2Theta(w);
}

inline float CosPhi(const float3& w) {
    float sinTheta = SinTheta(w);
    return (sinTheta == 0) ? 1 : clamp(w.x / sinTheta, -1.f, 1.f);
}

inline float SinPhi(const float3& w) {
    float sinTheta = SinTheta(w);
    return (sinTheta == 0) ? 0 : clamp(w.y / sinTheta, -1.f, 1.f);
}

inline float Cos2Phi(const float3& w) { return CosPhi(w) * CosPhi(w); }

inline float Sin2Phi(const float3& w) { return SinPhi(w) * SinPhi(w); }

static void TrowbridgeReitzSample11(float cosTheta, float U1, float U2, float* slope_x, float* slope_y) {
    // special case (normal incidence)
    if (cosTheta > .9999f) {
        float r = sqrt(U1 / (1 - U1));
        float phi = 6.28318530718f * U2;
        *slope_x = r * cos(phi);
        *slope_y = r * sin(phi);
        return;
    }

    float sinTheta = sqrt(max(0.f, 1.f - cosTheta * cosTheta));
    float tanTheta = sinTheta / cosTheta;
    float a = 1 / tanTheta;
    float G1 = 2 / (1 + std::sqrt(1.f + 1.f / (a * a)));

    // sample slope_x
    float A = 2 * U1 / G1 - 1;
    float tmp = 1.f / (A * A - 1.f);
    if (tmp > 1e10) tmp = 1e10;
    float B = tanTheta;
    float D = sqrt(max(float(B * B * tmp * tmp - (A * A - B * B) * tmp), 0.f));
    float slope_x_1 = B * tmp - D;
    float slope_x_2 = B * tmp + D;
    *slope_x = (A < 0 || slope_x_2 > 1.f / tanTheta) ? slope_x_1 : slope_x_2;

    // sample slope_y
    float S;
    if (U2 > 0.5f) {
        S = 1.f;
        U2 = 2.f * (U2 - .5f);
    }
    else {
        S = -1.f;
        U2 = 2.f * (.5f - U2);
    }
    float z =
        (U2 * (U2 * (U2 * 0.27385f - 0.73369f) + 0.46341f)) /
        (U2 * (U2 * (U2 * 0.093073f + 0.309420f) - 1.000000f) + 0.597999f);
    *slope_y = S * z * std::sqrt(1.f + *slope_x * *slope_x);
}

static float3 TrowbridgeReitzSample(const float3& wi, float alpha_x, float alpha_y, float U1, float U2) {
    // 1. stretch wi
    float3 wiStretched = normalize(float3(alpha_x * wi.x, alpha_y * wi.y, wi.z));

    // 2. simulate P22_{wi}(x_slope, y_slope, 1, 1)
    float slope_x, slope_y;
    TrowbridgeReitzSample11(CosTheta(wiStretched), U1, U2, &slope_x, &slope_y);

    // 3. rotate
    float tmp = CosPhi(wiStretched) * slope_x - SinPhi(wiStretched) * slope_y;
    slope_y = SinPhi(wiStretched) * slope_x + CosPhi(wiStretched) * slope_y;
    slope_x = tmp;

    // 4. unstretch
    slope_x = alpha_x * slope_x;
    slope_y = alpha_y * slope_y;

    // 5. compute normal
    return normalize(float3(-slope_x, -slope_y, 1.));
}

// TrowbridgeReitzDistribution
class MicrofacetDistribution {
public:
    static inline float RoughnessToAlpha(float roughness) {
        roughness = max(roughness, (float)1e-3);
        float x = log(roughness);
        return 1.62142f + 0.819955f * x + 0.1734f * x * x + 0.0171201f * x * x * x + 0.000640711f * x * x * x * x;
    }

    MicrofacetDistribution(float alphax, float alphay) : 
        alphax(std::max(0.001f, alphax)), 
        alphay(std::max(0.001f, alphay)) {}

    float D(const float3& wh) const {
        float tan2Theta = Tan2Theta(wh);
        if (isinf(tan2Theta)) return 0.;
        const float cos4Theta = Cos2Theta(wh) * Cos2Theta(wh);
        float e =
            (Cos2Phi(wh) / (alphax * alphax) + Sin2Phi(wh) / (alphay * alphay)) *
            tan2Theta;
        return 1 / (PI * alphax * alphay * cos4Theta * (1 + e) * (1 + e));
    }

    float3 Sample_wh(const float3& wo, const float2& u) const {
        float3 wh;
        bool flip = wo.z < 0;
        wh = TrowbridgeReitzSample(flip ? -wo : wo, alphax, alphay, u[0], u[1]);
        if (flip) wh = -wh;
        return wh;
    }

    float G1(const float3& w) const {
        return 1 / (1 + Lambda(w));
    }

    float G(const float3& wo, const float3& wi) const {
        return 1 / (1 + Lambda(wo) + Lambda(wi));
    }

    float Pdf(const float3& wo, const float3& wh) const {
        return D(wh) * G1(wo) * absdot(wo, wh) / AbsCosTheta(wo);
    }

private:
    float Lambda(const float3& w) const {
        float absTanTheta = abs(TanTheta(w));
        if (isinf(absTanTheta)) return 0.f;
        // Compute _alpha_ for direction _w_
        float alpha = sqrt(Cos2Phi(w) * alphax * alphax + Sin2Phi(w) * alphay * alphay);
        float alpha2Tan2Theta = (alpha * absTanTheta) * (alpha * absTanTheta);
        return (-1 + sqrt(1.f + alpha2Tan2Theta)) / 2;
    }

    const float alphax, alphay;
};

// https://seblagarde.wordpress.com/2013/04/29/memo-on-fresnel-equations/
float3 FrConductor(float cosThetaI, const float3& etai, const float3& etat, const float3& k) {
    cosThetaI = clamp(cosThetaI, -1.f, 1.f);
    float3 eta = etat / etai;
    float3 etak = k / etai;

    float cosThetaI2 = cosThetaI * cosThetaI;
    float sinThetaI2 = 1.f - cosThetaI2;
    float3 eta2 = eta * eta;
    float3 etak2 = etak * etak;

    float3 t0 = eta2 - etak2 - sinThetaI2;
    float3 a2plusb2 = sqrt(t0 * t0 + 4 * eta2 * etak2);
    float3 t1 = a2plusb2 + cosThetaI2;
    float3 a = sqrt(0.5f * (a2plusb2 + t0));
    float3 t2 = 2.f * cosThetaI * a;
    float3 Rs = (t1 - t2) / (t1 + t2);

    float3 t3 = cosThetaI2 * a2plusb2 + sinThetaI2 * sinThetaI2;
    float3 t4 = t2 * sinThetaI2;
    float3 Rp = Rs * (t3 - t4) / (t3 + t4);

    return 0.5 * (Rp + Rs);
}

float FrDielectric(float cosThetaI, float etaI, float etaT) {
    cosThetaI = clamp(cosThetaI, -1.f, 1.f);
    // Potentially swap indices of refraction
    bool entering = cosThetaI > 0.f;
    if (!entering) {
        swap(etaI, etaT);
        cosThetaI = abs(cosThetaI);
    }

    // Compute _cosThetaT_ using Snell's law
    float sinThetaI = sqrt(max(0.f, 1.f - cosThetaI * cosThetaI));
    float sinThetaT = etaI / etaT * sinThetaI;

    // Handle total internal reflection
    if (sinThetaT >= 1) return 1;
    float cosThetaT = sqrt(max(0.f, 1.f - sinThetaT * sinThetaT));
    float Rparl = ((etaT * cosThetaI) - (etaI * cosThetaT)) /
        ((etaT * cosThetaI) + (etaI * cosThetaT));
    float Rperp = ((etaI * cosThetaI) - (etaT * cosThetaT)) /
        ((etaI * cosThetaI) + (etaT * cosThetaT));
    return (Rparl * Rparl + Rperp * Rperp) / 2;
}

class Fresnel {
public:
    // Fresnel Interface
    virtual ~Fresnel() {}
    virtual float3 Evaluate(float cosI) const = 0;
};

class FresnelConductor : public Fresnel {
public:
    FresnelConductor(const float3& etaI, const float3& etaT, const float3& k) : etaI(etaI), etaT(etaT), k(k) {}
    virtual float3 Evaluate(float cosThetaI) const override {
        return FrConductor(abs(cosThetaI), etaI, etaT, k);
    }
private:
    float3 etaI, etaT, k;
};

class FresnelDielectric : public Fresnel {
public:
    FresnelDielectric(float etaI, float etaT) : etaI(etaI), etaT(etaT) {}
    virtual float3 Evaluate(float cosThetaI) const override {
        return FrDielectric(cosThetaI, etaI, etaT);
    }
private:
    float etaI, etaT;
};


class MicrofacetReflection {
public:
    MicrofacetReflection(const float3& R, MicrofacetDistribution* distribution, Fresnel* fresnel)
        : R(R), distribution(distribution), fresnel(fresnel) {}
    float3 f(const float3& wo, const float3& wi) const {
        float cosThetaO = AbsCosTheta(wo), cosThetaI = AbsCosTheta(wi);
        float3 wh = wi + wo;
        // Handle degenerate cases for microfacet reflection
        if (cosThetaI == 0 || cosThetaO == 0) return float3(0.f);
        if (wh.x == 0 && wh.y == 0 && wh.z == 0) return float3 (0.f);
        wh = normalize(wh);
        // For the Fresnel call, make sure that wh is in the same hemisphere
        // as the surface normal, so that TIR is handled correctly.
        float3 F = fresnel->Evaluate(dot(wi, Faceforward(wh, float3(0, 0, 1))));
        return R * distribution->D(wh) * distribution->G(wo, wi) * F /
            (4 * cosThetaI * cosThetaO);
    }
    float3 Sample_f(const float3& wo, float3* wi, const float2& u, float* pdf) const {
        // Sample microfacet orientation $\wh$ and reflected direction $\wi$
        if (wo.z == 0) return 0.;
        float3 wh = distribution->Sample_wh(wo, u);
        if (dot(wo, wh) < 0) return 0.;   // Should be rare
        *wi = reflect(wo, wh);
        if (!SameHemisphere(wo, *wi)) return float3(0.f);

        // Compute PDF of _wi_ for microfacet reflection
        *pdf = distribution->Pdf(wo, wh) / (4 * dot(wo, wh));
        return f(wo, *wi);
    }
    float Pdf(const float3& wo, const float3& wi) const {
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