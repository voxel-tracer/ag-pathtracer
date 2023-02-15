#pragma once

#include "rendertimer.h"

namespace Tmpl8
{

class Accumulator {
public:
	Accumulator(int w, int h) : width(w), height(h), samples(0) {
		pixels = (float3*)MALLOC64(width * height * sizeof(float3));
		Clear();
	}

	~Accumulator() { FREE64(pixels); }

	inline void AddSample(int x, int y, const float3& clr) {
		pixels[y * width + x] += clr;
	}

	inline void IncrementSampleCount() {
		samples++;
	}

	inline void Clear() {
		const int s = width * height;
		for (int i = 0; i < s; i++) pixels[i] = float3(0.f);

		samples = 0;
	}

	inline void CopyToSurface(Surface* screen, int offsetX, int offsetY) const {
		for (auto y = 0; y < height; y++) {
			for (auto x = 0; x < width; x++) {
				screen->Plot(x + offsetX, y + offsetY, rgb2uint(pixels[y * width + x] / (float)samples));
			}
		}
	}

private:
	float3* pixels;
	int width, height;
	int samples;
};


class MyApp : public TheApp
{
public:
	// game flow methods
	void Init();
	void Tick( float deltaTime );
	void Shutdown() { /* implement if you want to do something on exit */ }
	// input handling
	void MouseUp( int button ) { 
		if (button == GLFW_MOUSE_BUTTON_1)
			leftMouseBtnDown = false;
	}
	void MouseDown( int button ) { 
		if (button == GLFW_MOUSE_BUTTON_1)
			leftMouseBtnDown = true;
	}
	void MouseMove( int x, int y ) { 
		if (leftMouseBtnDown) {
			mouseDiff.x = x - mousePos.x;
			mouseDiff.y = y - mousePos.y;
		}
		mousePos.x = x, mousePos.y = y;
	}
	void MouseWheel( float y ) { /* implement if you want to handle the mouse wheel */ }
	void KeyUp( int key ) { /* implement if you want to handle keys */ }
	void KeyDown( int key ) { /* implement if you want to handle keys */ }
	// data members
	int2 mousePos;
	int2 mouseDiff;
	bool leftMouseBtnDown;

	shared_ptr<RotatingCamera> camera;
	shared_ptr<Scene> scene;
	shared_ptr<Integrator> integrator;

	RenderTimer timer;
	shared_ptr<Accumulator> accumulator;
};

} // namespace Tmpl8