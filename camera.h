#pragma once

class Ray {
public:
	Ray() = default;
	Ray(float3 o, float3 d, float t = FLT_MAX) : O(o), D(normalize(d)), t(t) {}

	float3 at(float t) const {
		return O + t * D;
	}
public:
	float3 O;
	float3 D;
	mutable float t;
};

struct CameraDesc {
	float3 lookfrom;		// camera position
	float3 lookat;			// point that will show up in the center of the screen
	float3 vup;				// used to compute up direction of camera coordinate system
	float aspect_ratio;		// screen width/size
	float vfov = 45;		// vertical field of view in degrees
	float focus_dist = 1.0f;	// distance of the screen from the camera origin (only matters once we support depth of field)
	float aperture = 0.0f;
};

class Camera {
public:
	Camera(const CameraDesc& desc) : Camera(
		desc.lookfrom,
		desc.lookat,
		desc.vup,
		desc.aspect_ratio,
		desc.vfov,
		desc.aperture
	) {}

	Camera(
		float3 lookfrom,
		float3 lookat,
		float3 vup,
		float aspect_ratio,
		float vfov,
		float aperture
	) : lookat(lookat), vup(vup) {
		auto theta = radians(vfov);

		// compute screen/viewport size assuming it's positioned 1 unit in from of the camera
		auto h = tan(theta / 2);
		viewport_height = 2 * h;
		viewport_width = aspect_ratio * viewport_height;

		lens_radius = aperture / 2;

		updateCoords(lookfrom);
	}

	Ray GetRay(float s, float t) const {
		float3 rd(0.f);
		if (lens_radius > 0.f) rd = lens_radius * RandomInUnitDisk();
		auto offset = u * rd.x + v * rd.y;
		float3 pixel = lower_left_corner + s * horizontal + t * vertical;
		return Ray(origin + offset, pixel - origin - offset);
	}

protected:
	void updateCoords(float3 lookfrom) {
		focus_dist = length(lookat - lookfrom);

		// compute camera coordinate system vectors
		w = normalize(lookfrom - lookat); // notice w points opposite to the screen
		u = normalize(cross(vup, w));
		v = cross(w, u);
		origin = lookfrom;

		// compute helper vectors used to compute pixels on the screen surface
		horizontal = focus_dist * viewport_width * u;
		vertical = focus_dist * viewport_height * v;
		lower_left_corner = origin - horizontal / 2 - vertical / 2 - focus_dist * w;
	}

protected:
	float3 lookat;
	float3 vup;
	float focus_dist;
	float lens_radius;

	float viewport_width;
	float viewport_height;

	float3 origin;
	float3 u, v, w;

	float3 lower_left_corner;
	float3 horizontal;
	float3 vertical;
};

class RotatingCamera : public Camera {
public:
	RotatingCamera(const CameraDesc& desc) : RotatingCamera(
		desc.lookfrom,
		desc.lookat,
		desc.vup,
		desc.aspect_ratio,
		desc.vfov,
		desc.aperture
	) {}
	RotatingCamera(
		float3 lookfrom,
		float3 lookat,
		float3 vup,
		float aspect_ratio,
		float vfov,
		float aperture
	) : Camera(
		lookfrom,
		lookat,
		vup,
		aspect_ratio,
		vfov,
		aperture
	) {
		dist2LookAt = length(lookat - lookfrom);
		// compute rotations
		float3 lfn = normalize(lookfrom - lookat);
		float xa = acosf(sqrtf(lfn.x * lfn.x + lfn.z * lfn.z));
		if (lfn.y > 0) xa = -xa;
		float ya = acosf(lfn.z / sqrtf(lfn.x * lfn.x + lfn.z * lfn.z));
		if (lfn.x < 0) ya = -ya;

		xAngle = xa;
		yAngle = ya;
	}

	void update(float2 angle) {
		xAngle += angle.x;
		xAngle = clamp(xAngle, 0.0f, PI / 2 - 0.01f);

		yAngle += angle.y;

		mat4 mat = mat4::RotateY(yAngle) * mat4::RotateX(xAngle);

		float3 lookfrom = mat.TransformVector(float3(0, 0, 1)) * dist2LookAt + lookat;
		updateCoords(lookfrom);
	}
private:
	//float3 lookfrom;
	float dist2LookAt;
	float xAngle;
	float yAngle;
};
