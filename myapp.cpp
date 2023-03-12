#include "precomp.h"
#include "disney.h"
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

	auto mat = Material::make_disney(rgb2lin(float3(.529f, .145f, .039f)), .25f, 0.f);
	auto floor = Material::make_lambertian(checker);

	vector<shared_ptr<Intersectable>> primitives;
	primitives.push_back(std::make_shared<Plane>(make_float3(0, -1, 0), make_float2(20), floor));
	mat4 transform = mat4::Translate(0, -1, 0) * mat4::RotateY(radians(60));
	
	auto trimesh = TriangleMesh::LoadObj("D://models/bunny.obj", mat, transform);
	primitives.push_back(make_shared<BVHTriMesh>(trimesh, 1));

	CameraDesc camera;
	camera.lookfrom = float3(3, 1.5f, 4);
	camera.lookat = float3(.5f, 0, .5f);
	camera.vup = float3(0, 1, 0);
	camera.aspect_ratio = 1;
	camera.vfov = 30;
	//camera.aperture = .1f;

	// add an emitting sphere
	auto lightE = float3(523, 342, 342);
	auto lightMat = Material::make_emitter(lightE);
	auto light = make_shared<Sphere>(float3(-30, 100, 40), 5.f, lightMat);
	primitives.push_back(light);

	auto scene = make_shared<Scene>(primitives, camera);
	scene->lights.push_back(make_shared<AreaLight>(light, lightE));
	scene->lights.push_back(make_shared<InfiniteAreaLight>(float3(.4f, .45f, .5f)));

	return scene;
}

shared_ptr<Intersectable> makeSphere(const float3& center, float radius, const float3 color, float roughness) {
	shared_ptr<Material> mat = Material::make_disney(color, roughness, 0.f);
	return make_shared<Sphere>(center, radius, mat);
}

// Keep this as it matches PBRT's simple.pbrt
shared_ptr<Scene> SimpleTestScene() {
	auto light_grey = make_shared<SolidColor>(.8f);
	auto dark_grey = make_shared<SolidColor>(.1f);
	auto checker = make_shared<CheckerTexture>(light_grey, dark_grey);
	auto floor = Material::make_lambertian(checker);

	auto disneyDielectric = Material::make_disney(rgb2lin(float3(.529f, .145f, .039f)), .25f, 0.f);

	vector<shared_ptr<Intersectable>> primitives;
	primitives.push_back(std::make_shared<Plane>(make_float3(0, -1, 0), make_float2(20), floor));
	//primitives.push_back(make_shared<Sphere>(float3(0.f), 1.f, disneyDielectric));

	// add an emitting sphere
	auto lightE = float3(523, 342, 342);
	auto lightMat = Material::make_emitter(lightE);
	auto light = make_shared<Sphere>(float3(-30, 100, 40), 5.f, lightMat);
	primitives.push_back(light);

	CameraDesc camera;
	camera.lookfrom = float3(3.f, 1.5f, 4.f);
	camera.lookat = float3(0, 0, 0);
	camera.vup = float3(0, 1, 0);
	camera.aspect_ratio = 1;
	//camera.aperture = .1f;

	auto scene = make_shared<Scene>(primitives, camera);
	scene->lights.push_back(make_shared<InfiniteAreaLight>(float3(.4f, .45f, .5f)));
	scene->lights.push_back(make_shared<AreaLight>(light, lightE));

	return scene;
}

TheApp* CreateApp() { return new MyApp(); }

// -----------------------------------------------------------
// Initialize the application
// -----------------------------------------------------------
void MyApp::Init()
{
	// anything that happens only once at application start goes here
	scene = SimpleTestScene();
	camera = make_shared<RotatingCamera>(scene->camera);

	integratorR = make_shared<PathTracer>();
	integratorL = integratorR;

	const int MinScrSize = min(SCRWIDTH, SCRHEIGHT);
	int2 scrPos((SCRWIDTH - MinScrSize) / 2, (SCRHEIGHT - MinScrSize) / 2);
	accumulator = make_shared<Accumulator>(MinScrSize, MinScrSize, scrPos);
}

// -----------------------------------------------------------
// Main application tick function - Executed once per frame
// deltaTime: elapsed time since last frame in milliseconds
// -----------------------------------------------------------
void MyApp::Tick( float deltaTime )
{
	// rotate camera if mouse moved since last tick
	if (!paused && (mouseDiff.x != 0 || mouseDiff.y != 0)) {
		float rotation_speed = 0.005f;
		camera->update(float2(-mouseDiff.y * rotation_speed, -mouseDiff.x * rotation_speed));
		accumulator->Clear();
	}
	mouseDiff = int2(0); // reset mouse diff

	// clear the screen to black
	screen->Clear(0);

	if (!paused) {
		timer.start();

		accumulator->IncrementSampleCount();

		for (int y = 0; y < accumulator->height; y++) {
			for (int x = 0; x < accumulator->width; x++) {
				float2 p(x + RandomFloat(), y + RandomFloat());
				float2 uv = accumulator->PixelToFilm(p);
				Ray ray = camera->GetRay(uv.x, uv.y);
				float3 clr = (uv.x < 0.5f) ? integratorL->Li(ray, *scene) : integratorR->Li(ray, *scene, 0, true);
				if (HasNans(clr) || std::isinf(Luminance(clr))) {
					//std::cerr << "outlier\n";
					clr = float3(0.f);
				}
				accumulator->AddSample(x, y, clr);
			}
		}

		if (timer.stop()) {
			auto stats = timer.getStatsAndReset();
			printf("Frame %d\tRender stats { best = %.2f, avg = %.2f, worst = %.2f }\t\t\r",
				accumulator->NumSamples(), stats.bestDuration, stats.avgDuration, stats.worstDuration);
		}
	}


	accumulator->CopyToSurface(screen);

	if (paused && accumulator->IsInside(mousePos.x, mousePos.y)) {
		auto fc = accumulator->WindowToFilm(mousePos);
		Ray ray = camera->GetRay(fc.x, fc.y);
		SurfaceInteraction hit;
		if (scene->NearestIntersection(ray, hit)) {
			auto p = accumulator->FilmToWindow(camera->WorldToScreen(hit.p));
			auto n = accumulator->FilmToWindow(camera->WorldToScreen(hit.p + hit.shading.n));
			auto dpdu = accumulator->FilmToWindow(camera->WorldToScreen(hit.p + hit.dpdu));
			auto dpdv = accumulator->FilmToWindow(camera->WorldToScreen(hit.p + hit.dpdv));
			screen->Line(p.x, p.y, n.x, n.y, 0x0000FF);
			screen->Line(p.x, p.y, dpdu.x, dpdu.y, 0xFF0000);
			screen->Line(p.x, p.y, dpdv.x, dpdv.y, 0x00FF00);
			screen->Box((int)p.x - 5, (int)p.y - 5, (int)p.x + 5, (int)p.y + 5, 0xFF0000);
		}
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