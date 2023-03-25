#include "precomp.h"

#include "intersectable.h"
#include "lights.h"

float3 SurfaceInteraction::Le(const float3& w) const {
	const AreaLight* area = shape->GetAreaLight();
	return area ? area->L : float3(0.f);
}