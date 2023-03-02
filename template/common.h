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


// IMPORTANT NOTE ON OPENCL COMPATIBILITY ON OLDER LAPTOPS:
// Without a GPU, a laptop needs at least a 'Broadwell' Intel CPU (5th gen, 2015):
// Intel's OpenCL implementation 'NEO' is not available on older devices.
// Same is true for Vulkan, OpenGL 4.0 and beyond, as well as DX11 and DX12.