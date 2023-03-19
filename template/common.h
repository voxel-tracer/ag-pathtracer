// common.h is to be included in host and device code and stores
// global settings and defines.

// default screen resolution
#define SCRWIDTH	400
#define SCRHEIGHT	400
// #define FULLSCREEN

// constants
#define PI		3.14159265358979323846264f
#define INVPI		0.31830988618379067153777f
#define INV2PI		0.15915494309189533576888f
#define TWOPI		6.28318530717958647692528f
#define SQRT_PI_INV	0.56418958355f
#define LARGE_FLOAT	1e34f
#define EPSILON		0.0001f

static const float FloatOneMinusEpsilon = 0x1.fffffep-1;

// Utility functions
inline float radians(float degrees) {
	return degrees * PI / 180.0f;
}

inline float degrees(float radians) {
	return radians * 180.0f * INVPI;
}

inline float3 rgb2lin(float3 c) {
	return { pow(c.x, 2.2f), pow(c.y, 2.2f),pow(c.z, 2.2f) };
}

inline float3 hex2lin(int hexValue) {
	float3 rgb;
	rgb.x = ((hexValue >> 16) & 0xFF) / 255.f;
	rgb.y = ((hexValue >> 8) & 0xFF) / 255.f;
	rgb.z = (hexValue & 0xFF) / 255.f;
	return rgb2lin(rgb);
}

inline float3 lin2rgb(float3 c) {
	float e = 1 / 2.2f;
	return { pow(c.x, e), pow(c.y, e),pow(c.z, e) };
}

inline uint rgb2uint(float3 clr) {
	int r = static_cast<int>(256 * clamp(clr.x, 0.0, 0.999));
	int g = static_cast<int>(256 * clamp(clr.y, 0.0, 0.999));
	int b = static_cast<int>(256 * clamp(clr.z, 0.0, 0.999));
	return (r << 16) + (g << 8) + b;
}

inline float GammaCorrect(float value) {
	if (value <= 0.0031308f) return 12.92f * value;
	return 1.055f * std::pow(value, (1.f / 2.4f)) - 0.055f;
}


inline float RandomFloat(float min, float max) {
	// Returns a random real in [min,max).
	return min + (max - min) * RandomFloat();

}

inline float3 RandomInUnitDisk() {
	while (true) {
		auto p = float3(RandomFloat(-1, 1), RandomFloat(-1, 1), 0);
		if (sqrLength(p) >= 1) continue;
		return p;
	}
}

inline float3 RandomInSphere(float Radius = 1.f) {
	auto a = 1 - 2 * RandomFloat();
	auto b = sqrt(1 - a * a);
	auto phi = 2 * PI * RandomFloat();
	return make_float3(
		Radius * b * cos(phi),
		Radius * b * sin(phi),
		Radius * a
	);
}

inline float3 RandomInSphere(const float2& u) {
	auto a = 1 - 2 * u[0];
	auto b = sqrt(1 - a * a);
	auto phi = 2 * PI * u[1];
	return make_float3(b * cos(phi), b * sin(phi), a);
}

inline float3 RandomInHemisphere(const float3& N) {
	auto in_unit_sphere = RandomInSphere();
	// if point is not in same hemisphere as normal N, invert it
	if (dot(in_unit_sphere, N) < 0)
		return -in_unit_sphere;
	return in_unit_sphere;
}

inline float3 TangentToWorld(const float3& N, const float3& V) {
	auto w = N;
	auto a = (fabs(w.x > .9f)) ? float3(0, 1, 0) : float3(1, 0, 0);
	auto v = normalize(cross(w, a));
	auto u = cross(w, v);
	return V.x * u + V.y * v + V.z * w;
}

inline float3 CosineWeightedRandomInHemisphere(const float3& N) {
	// work in tangent space where N=(0, 0, 1)
	float r0 = RandomFloat(), r1 = RandomFloat();
	float r = sqrt(r0);
	float theta = TWOPI * r1;
	float x = r * cos(theta);
	float y = r * sin(theta);
	float3 R(x, y, sqrt(1 - r0));
	return TangentToWorld(N, R);
}

inline float2 ConcentricSampleDisk(const float2& u) {
	// Map uniform random numbers to [-1,1]^2
	float2 uOffset = 2.f * u - float2(1, 1);

	// Handle degeneracy at the origin
	if (uOffset.x == 0 && uOffset.y == 0) return float2(0, 0);

	// Apply concentric mapping to point
	float theta, r;
	if (std::abs(uOffset.x) > std::abs(uOffset.y)) {
		r = uOffset.x;
		theta = (PI / 4) * (uOffset.y / uOffset.x);
	}
	else {
		r = uOffset.y;
		theta = (PI / 2) - (PI / 4) * (uOffset.x / uOffset.y);
	}
	return r * float2(std::cos(theta), std::sin(theta));
}

// Cosine sampling of hemisphere in shading space (normal is (0, 0, 1))
inline float3 CosineSampleHemisphere(const float2& u) {
	float2 d = ConcentricSampleDisk(u);
	float z = std::sqrt(std::max(0.f, 1 - d.x * d.x - d.y * d.y));
	return float3(d.x, d.y, z);
}

inline void CoordinateSystem(const float3& v1, float3* v2, float3* v3) {
	if (std::abs(v1.x) > std::abs(v1.y))
		*v2 = float3(-v1.z, 0, v1.x) / std::sqrt(v1.x * v1.x + v1.z * v1.z);
	else
		*v2 = float3(0, v1.z, -v1.y) / std::sqrt(v1.y * v1.y + v1.z * v1.z);
	*v3 = cross(v1, *v2);
}

inline float3 SphericalDirection(float sinTheta, float cosTheta, float phi,
		const float3& x, const float3& y, const float3& z) {
	return sinTheta * std::cos(phi) * x + sinTheta * std::sin(phi) * y + cosTheta * z;
}

inline float UniformConePdf(float cosThetaMax) {
	return 1 / (2 * PI * (1 - cosThetaMax));
}

// IMPORTANT NOTE ON OPENCL COMPATIBILITY ON OLDER LAPTOPS:
// Without a GPU, a laptop needs at least a 'Broadwell' Intel CPU (5th gen, 2015):
// Intel's OpenCL implementation 'NEO' is not available on older devices.
// Same is true for Vulkan, OpenGL 4.0 and beyond, as well as DX11 and DX12.