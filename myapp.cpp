#include "precomp.h"
#include "whitted.h"
#include "bvhtrimesh.h"
#include "myapp.h"

#define TINYOBJLOADER_IMPLEMENTATION
#define TINYOBJLOADER_USE_MAPBOX_EARCUT
#include "tiny_obj_loader.h"

shared_ptr<Scene> BunnyScene() {
	float3 dark_red(.1f, 0.f, 0.f);
	auto light_grey = make_shared<SolidColor>(.8f);
	auto dark_grey = make_shared<SolidColor>(.1f);
	auto checker = make_shared<CheckerTexture>(light_grey, dark_grey);

	auto mat = Material::make_lambertian(SolidColor::make(dark_red));
	//auto mat = Material::make_glass(1.05f);
	auto floor = Material::make_lambertian(checker);

	vector<shared_ptr<Intersectable>> primitives;
	primitives.push_back(std::make_shared<Plane>(make_float3(0, 1, 0), make_float2(20), floor));
	mat4 transform = mat4::Translate(0, 1, 0) * mat4::RotateY(radians(180)) * mat4::RotateX(radians(180));
	
	auto trimesh = TriangleMesh::LoadObj("D://models/bunny-0.1.obj", mat, transform, true);
	primitives.push_back(make_shared<BVHTriMesh>(trimesh, 1));
	//primitives.push_back(trimesh);

	CameraDesc camera{ { 3.f, -1.5f, 4.f }, { .5f, 0, .5f }, { 0.f, 1.f, 0.f }, 1.f };
	float3 light = float3(-30, -100, 40);
	float3 lightColor = float3(52300, 34200, 34200); // approximately 4000K black body light source

	auto scene = make_shared<Scene>(primitives, camera);
	scene->lights.push_back(make_shared<InfiniteAreaLight>(float3(.4f, .45f, .5f)));
	scene->lights.push_back(make_shared<PointLight>(light, lightColor));

	return scene;
}

shared_ptr<Scene> TestGlassScene() {
	float3 dark_red(.1f, 0.f, 0.f);
	auto dark_grey = make_shared<SolidColor>(.8f);
	auto light_grey = make_shared<SolidColor>(.1f);
	auto checker = make_shared<CheckerTexture>(light_grey, dark_grey);

	auto glass = Material::make_glass(1.05f, dark_red);
	auto floor = Material::make_lambertian(checker);

	vector<shared_ptr<Intersectable>> primitives;
	primitives.push_back(std::make_shared<Plane>(make_float3(0, 1, 0), make_float2(20), floor));
	primitives.push_back(std::make_shared<Sphere>(make_float3(-2.0f, 0, 0), 1.0f, glass));
	primitives.push_back(std::make_shared<Sphere>(make_float3(0, 0.5f, 0), 0.5f, glass));
	primitives.push_back(std::make_shared<Sphere>(make_float3(1.0f, 0.75f, 0), 0.25f, glass));

	CameraDesc camera{ { 3.f, -1.5f, 4.f }, { .5f, 0, .5f }, { 0.f, 1.f, 0.f }, 1.f };
	float3 light = float3(-30, -100, 40);
	float3 lightColor = float3(52300, 34200, 34200); // approximately 4000K black body light source

	auto scene = make_shared<Scene>(primitives, camera);
	scene->lights.push_back(make_shared<InfiniteAreaLight>(float3(.4f, .45f, .5f)));

	return scene;
}

//Scene InitScene() {
//	auto orange = make_shared<SolidColor>(rgb2lin({ 1.f, .68f, .25f }));
//	auto yellow = make_shared<SolidColor>(rgb2lin({ .98f, .85f, .37f }));
//	auto checker = make_shared<CheckerTexture>(yellow, orange);
//	auto grey = make_shared<SolidColor>(.73f);
//
//	auto glass = Material::make_glass(1.05f);
//	auto white = make_shared<Material>(grey);
//	auto mirror = Material::make_mirror(.8f, .8f, .8f);
//	auto floor = make_shared<Material>(checker);
//
//	vector<shared_ptr<Intersectable>> primitives;
//	primitives.push_back(std::make_shared<Plane>(make_float3(0, 1, 0), make_float2(4, 4), floor));
//	primitives.push_back(std::make_shared<Sphere>(make_float3(.5f, .2f, -.25f), .5f, glass));
//	primitives.push_back(std::make_shared<Sphere>(make_float3(-.65f, 0, -.75f), .45f, mirror));
//
//	CameraDesc camera{ { .25f, -1, 2.5 }, { .25f, 0.f, 0.f }, { 0.f, 1.f, 0.f }, 1.f, 45 };
//	float3 light{ 0.0f, -10, 0.0f };
//	float3 lightColor(100);
//
//	return Scene(primitives, light, lightColor, camera);
//}

TheApp* CreateApp() { return new MyApp(); }

// -----------------------------------------------------------
// Initialize the application
// -----------------------------------------------------------
void MyApp::Init()
{
	// anything that happens only once at application start goes here
	//InitScene();
	//scene = TestGlassScene();
	scene = BunnyScene();
	camera = make_shared<RotatingCamera>(scene->camera);
	integrator = make_shared<WhittedIntegrator>(10);
}

// -----------------------------------------------------------
// Main application tick function - Executed once per frame
// deltaTime: elapsed time since last frame in milliseconds
// -----------------------------------------------------------
void MyApp::Tick( float deltaTime )
{
	// rotate camera if mouse moved since last tick
	if (mouseDiff.x != 0 || mouseDiff.y != 0) {
		float rotation_speed = 0.005f;
		camera->update(float2(mouseDiff.y * rotation_speed, -mouseDiff.x * rotation_speed));
	}
	mouseDiff = int2(0); // reset mouse diff

	// clear the screen to black
	screen->Clear(0);

	const int MinScrSize = min(SCRWIDTH, SCRHEIGHT);

	timer.start();

	const int RenderWidth = MinScrSize;
	const int RenderHeight = MinScrSize;
	for (int y = -RenderHeight / 2; y < RenderHeight / 2; y++) {
		float v = ((float)y) / RenderHeight + .5f;
		for (int x = -RenderWidth / 2; x < RenderWidth / 2; x++) {
			float u = ((float)x) / RenderWidth + .5f;
			Ray ray = camera->GetRay(u, v);
			float3 clr = integrator->Li(ray, *scene);
			screen->Plot(SCRWIDTH / 2 + x, SCRHEIGHT / 2 + y, rgb2uint(clr));
		}
	}

	if (timer.stop()) {
		auto stats = timer.getStatsAndReset();
		printf("Render stats { best = %.2f, avg = %.2f, worst = %.2f }\t\t\r", stats.bestDuration, stats.avgDuration, stats.worstDuration);
	}

#if 0

	static Kernel* kernel = 0;			// statics should be members of MyApp of course.
	static Surface bitmap( 512, 512 );	// having them here allows us to disable the OpenCL
	static Buffer* clBuffer = 0;		// demonstration using a single #if 0.
	static int offset = 0;
	if (!kernel)
	{
		// prepare for OpenCL work
		Kernel::InitCL();		
		// compile and load kernel "render" from file "kernels.cl"
		kernel = new Kernel( "cl/kernels.cl", "render" );
		// create an OpenCL buffer over using bitmap.pixels
		clBuffer = new Buffer( 512 * 512, Buffer::DEFAULT, bitmap.pixels );
	}
	// pass arguments to the OpenCL kernel
	kernel->SetArgument( 0, clBuffer );
	kernel->SetArgument( 1, offset++ );
	// run the kernel; use 512 * 512 threads
	kernel->Run( 512 * 512 );
	// get the results back from GPU to CPU (and thus: into bitmap.pixels)
	clBuffer->CopyFromDevice();
	// show the result on screen
	bitmap.CopyTo( screen, 500, 200 );

#endif

}