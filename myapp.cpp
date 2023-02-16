#include "precomp.h"
#include "integrator.h"
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

	auto mat = Material::make_lambertian(SolidColor::make(float3(.7f, .1f, .1f)));
	//auto mat = Material::make_glass(1.05f);
	auto floor = Material::make_lambertian(checker);

	vector<shared_ptr<Intersectable>> primitives;
	primitives.push_back(std::make_shared<Plane>(make_float3(0, 1, 0), make_float2(20), floor));
	mat4 transform = mat4::Translate(0, 1, 0) * mat4::RotateY(radians(180)) * mat4::RotateX(radians(180));
	
	auto trimesh = TriangleMesh::LoadObj("D://models/bunny-lowpoly.obj", mat, transform, true);
	primitives.push_back(make_shared<BVHTriMesh>(trimesh, 1));

	CameraDesc camera{ { 3.f, -1.5f, 4.f }, { .5f, 0, .5f }, { 0.f, 1.f, 0.f }, 1.f };
	//camera.aperture = 1.f;

	// add an emitting sphere
	auto lightE = float3(523, 342, 342);
	auto lightMat = Material::make_emitter(523, 342, 342);
	auto light = make_shared<Sphere>(float3(-30, -100, 40), 10.f, lightMat);
	primitives.push_back(light);

	auto scene = make_shared<Scene>(primitives, camera);
	scene->lights.push_back(make_shared<AreaLight>(light, lightE));
	scene->lights.push_back(make_shared<InfiniteAreaLight>(float3(.4f, .45f, .5f)));
	//scene->lights.push_back(make_shared<PointLight>(float3(-30, -100, 40), float3(52300, 34200, 34200)));

	return scene;
}

shared_ptr<Scene> GlassScene() {
	float3 dark_red(.1f, 0.f, 0.f);
	auto light_grey = make_shared<SolidColor>(.8f);
	auto dark_grey = make_shared<SolidColor>(.1f);
	auto checker = make_shared<CheckerTexture>(light_grey, dark_grey);

	auto mat = Material::make_glass(1.125f);
	auto floor = Material::make_lambertian(checker);

	vector<shared_ptr<Intersectable>> primitives;
	primitives.push_back(std::make_shared<Plane>(make_float3(0, 1, 0), make_float2(20), floor));
	primitives.push_back(make_shared<Sphere>(float3(.5f, .25f - EPSILON, .5f), .75f, mat));

	CameraDesc camera;
	camera.lookfrom = float3(3.f, -1.5f, 4.f);
	camera.lookat = float3(.5f, 0, .5f);
	camera.vup = float3(0, 1, 0);
	camera.aspect_ratio = 1;
	//camera.aperture = .1f;

	auto scene = make_shared<Scene>(primitives, camera);
	scene->lights.push_back(make_shared<InfiniteAreaLight>(float3(.4f, .45f, .5f)));
	scene->lights.push_back(make_shared<PointLight>(float3(-30, -100, 40), float3(52300, 34200, 34200)));

	return scene;
}

TheApp* CreateApp() { return new MyApp(); }

// -----------------------------------------------------------
// Initialize the application
// -----------------------------------------------------------
void MyApp::Init()
{
	// anything that happens only once at application start goes here
	//scene = GlassScene();
	scene = BunnyScene();
	camera = make_shared<RotatingCamera>(scene->camera);

	//integrator = make_shared<WhittedIntegrator>(10);
	integratorL = make_shared<SimplePT2>();
	integratorR = make_shared<SimplePT>();

	const int MinScrSize = min(SCRWIDTH, SCRHEIGHT);
	accumulator = make_shared<Accumulator>(MinScrSize, MinScrSize);
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
		accumulator->Clear();
	}
	mouseDiff = int2(0); // reset mouse diff

	// clear the screen to black
	screen->Clear(0);

	const int MinScrSize = min(SCRWIDTH, SCRHEIGHT);

	timer.start();

	accumulator->IncrementSampleCount();

	const int RenderWidth = MinScrSize;
	const int RenderHeight = MinScrSize;
	for (int y = 0; y < RenderHeight; y++) {
		float v = (y + RandomFloat()) / RenderHeight;
		for (int x = 0; x < RenderWidth; x++) {
			float u = (x + RandomFloat()) / RenderWidth;
			Ray ray = camera->GetRay(u, v);
			float3 clr = (x < RenderWidth / 2) ? integratorL->Li(ray, *scene) : integratorR->Li(ray, *scene);
			accumulator->AddSample(x, y, clr);
		}
	}

	if (timer.stop()) {
		auto stats = timer.getStatsAndReset();
		printf("Frame %d\tRender stats { best = %.2f, avg = %.2f, worst = %.2f }\t\t\r", 
			accumulator->NumSamples(), stats.bestDuration, stats.avgDuration, stats.worstDuration);
	}

	accumulator->CopyToSurface(screen, (SCRWIDTH - RenderWidth) / 2, (SCRHEIGHT - RenderHeight) / 2);

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