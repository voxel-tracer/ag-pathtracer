// common.h is to be included in host and device code and stores
// global settings and defines.

// default screen resolution
#define SCRWIDTH	1280
#define SCRHEIGHT	720
// #define FULLSCREEN

// constants
#define PI		3.14159265358979323846264f
#define INVPI		0.31830988618379067153777f
#define INV2PI		0.15915494309189533576888f
#define TWOPI		6.28318530717958647692528f
#define SQRT_PI_INV	0.56418958355f
#define LARGE_FLOAT	1e34f
#define EPSILON		0.0001f

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


// IMPORTANT NOTE ON OPENCL COMPATIBILITY ON OLDER LAPTOPS:
// Without a GPU, a laptop needs at least a 'Broadwell' Intel CPU (5th gen, 2015):
// Intel's OpenCL implementation 'NEO' is not available on older devices.
// Same is true for Vulkan, OpenGL 4.0 and beyond, as well as DX11 and DX12.