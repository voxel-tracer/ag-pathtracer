#pragma once

#include "rendertimer.h"

namespace Tmpl8
{


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

	RenderTimer timer;
};

} // namespace Tmpl8