#include "precomp.h"
#include "disney.h"
#include "integrator.h"
#include "bvhtrimesh.h"
#include "texture.h"
#include "myapp.h"
#include "scene.h"

#define TINYOBJLOADER_IMPLEMENTATION
#define TINYOBJLOADER_USE_MAPBOX_EARCUT
#include "tiny_obj_loader.h"

void BunnyScene(Scene *scene) {
	auto red = DisneyMaterial::Make(rgb2lin(float3(.529f, .145f, .039f)), .25f, 0.f);
	auto gold = DisneyMaterial::Make((float3(0.944f, 0.776f, 0.373f)), 0.2f, 1.f);
	auto cute = DisneyMaterial::Make(hex2lin(0xc5b5d2), .25f, 0.f);

	auto floor = DisneyMaterial::Make(hex2lin(0xcbceb1), 1.f, 0.f);
	auto backdrop = TriangleMesh::CreateBackdrop(make_float3(0, -1, 20), float3(40, 20, 40), 7.5, 32, floor);
	scene->primitives.push_back(make_shared<BVHTriMesh>(backdrop, floor, 1));

	mat4 transform = mat4::Translate(.25, -1.05, .5) * mat4::RotateY(radians(180));
	auto trimesh = TriangleMesh::LoadObj("D://models/bunny.obj", cute, transform);
	scene->primitives.push_back(make_shared<BVHTriMesh>(trimesh, cute, 1));

	scene->camera.lookfrom = normalize(float3(1.24, 0.1, -4.84)) * 5;
	scene->camera.lookat = float3(.5f, 0, .5f);
	scene->camera.vup = float3(0, 1, 0);
	scene->camera.aspect_ratio = 1;
	scene->camera.vfov = 30;
	//scene->camera.aperture = .1f;


	// add area lights
	// key light
	scene->addAreaLight(
		make_shared<Sphere>(float3(0, 25, -20), 1.f, nullptr),
		float3(1, .941, .914) * 200
	);
	// fill light
	scene->addAreaLight(
		make_shared<Sphere>(float3(10, 25, -20), 1.f, nullptr),
		float3(1, .941, .914) * 50
	);
	// back light
	scene->addAreaLight(
		make_shared<Sphere>(float3(0, 20, 10), 5.f, nullptr),
		float3(1, .941, .914) * 1
	);

	scene->lights.push_back(make_shared<InfiniteAreaLight>(float3(.4f, .45f, .5f)));
}

// Keep this as it matches PBRT's simple.pbrt
void SimpleTestScene(Scene *scene) {

	//auto c = rgb2lin(float3(.529f, .145f, .039f));
	//std::cerr << c.x << ", " << c.y << ", " << c.z << std::endl;
	auto disneyDielectric = DisneyMaterial::Make(rgb2lin(float3(.529f, .145f, .039f)), .25f, 0.f);
	auto gold = DisneyMaterial::Make((float3(0.944f, 0.776f, 0.373f)), 0.25f, 1.f);
	auto aluminum = DisneyMaterial::Make((float3(0.912, 0.914, 0.920)), 0.25f, 0.f);
	auto bone = DisneyMaterial::Make((float3(0.793, 0.793, 0.664)), 0.25f, 0.f);
	auto brass = DisneyMaterial::Make((float3(0.887, 0.789, 0.434)), 0.25f, 0.f);
	auto brick = DisneyMaterial::Make((float3(0.262, 0.095, 0.061)), 0.25f, 0.f);
	auto charcoal = DisneyMaterial::Make((float3(0.020, 0.020, 0.020)), 0.25f, 0.f);
	auto chocolate = DisneyMaterial::Make((float3(0.162, 0.091, 0.060)), 0.25f, 0.f);
	auto chromium = DisneyMaterial::Make((float3(0.550, 0.556, 0.554)), 0.25f, 0.f);

	auto cute0 = DisneyMaterial::Make(hex2lin(0xf19a91), .25f, 0.f);
	auto cute1 = DisneyMaterial::Make(hex2lin(0xedd0ca), .25f, 0.f);
	auto cute2 = DisneyMaterial::Make(hex2lin(0xf3b8a8), .25f, 0.f);
	auto cute3 = DisneyMaterial::Make(hex2lin(0xf9ece6), .25f, 0.f);
	auto cute4 = DisneyMaterial::Make(hex2lin(0xf6e7d0), .25f, 0.f);
	auto cute5 = DisneyMaterial::Make(hex2lin(0xf5deac), .25f, 0.f);
	auto cute6 = DisneyMaterial::Make(hex2lin(0xeecf74), .25f, 0.f);
	auto cute7 = DisneyMaterial::Make(hex2lin(0x9ed5d8), .25f, 0.f);
	auto cute8 = DisneyMaterial::Make(hex2lin(0x9ba6ac), .25f, 0.f);
	auto cute9 = DisneyMaterial::Make(hex2lin(0xaebdc4), .25f, 0.f);
	auto cute10 = DisneyMaterial::Make(hex2lin(0xb9ddf3), .25f, 0.f);
	auto cute11 = DisneyMaterial::Make(hex2lin(0x87abc5), .25f, 0.f);
	auto cute12 = DisneyMaterial::Make(hex2lin(0xcbceb1), .25f, 0.f);
	auto cute13 = DisneyMaterial::Make(hex2lin(0xf7f7f7), .25f, 0.f);
	auto cute14 = DisneyMaterial::Make(hex2lin(0xc4ac64), .25f, 0.f);
	auto cute15 = DisneyMaterial::Make(hex2lin(0xe2f4f6), .25f, 0.f);
	auto cute16 = DisneyMaterial::Make(hex2lin(0xd2e4e6), .25f, 0.f);
	auto cute17 = DisneyMaterial::Make(hex2lin(0xbfdcda), .25f, 0.f);
	auto cute18 = DisneyMaterial::Make(hex2lin(0x69bab3), .25f, 0.f);
	auto cute19 = DisneyMaterial::Make(hex2lin(0x88cabc), .25f, 0.f);
	auto cute20 = DisneyMaterial::Make(hex2lin(0xcdd1d4), .25f, 0.f);
	auto cute21 = DisneyMaterial::Make(hex2lin(0xe6e5ea), .25f, 0.f);
	auto cute22 = DisneyMaterial::Make(hex2lin(0x33455b), .25f, 0.f);
	auto cute23 = DisneyMaterial::Make(hex2lin(0x5b6268), .25f, 0.f);
	auto cute24 = DisneyMaterial::Make(hex2lin(0x778592), .25f, 0.f);
	auto cute25 = DisneyMaterial::Make(hex2lin(0xe57a82), .25f, 0.f);
	auto cute26 = DisneyMaterial::Make(hex2lin(0xcd7d88), .25f, 0.f);
	auto cute27 = DisneyMaterial::Make(hex2lin(0xe3a3b1), .25f, 0.f);
	auto cute28 = DisneyMaterial::Make(hex2lin(0xf0d1e3), .25f, 0.f);
	auto cute = DisneyMaterial::Make(hex2lin(0xc5b5d2), .25f, 0.f);

	auto floor = DisneyMaterial::Make(hex2lin(0xcbceb1), 1.f, 0.f);
	auto backdrop = TriangleMesh::CreateBackdrop(make_float3(0, -1, 20), float3(40, 20, 40), 7.5, 32, floor);
	scene->primitives.push_back(make_shared<BVHTriMesh>(backdrop, floor, 1));

	scene->primitives.push_back(make_shared<Sphere>(float3(0, 0, 0), 1.f, gold));

	// add area lights
	// key light
	scene->addAreaLight(
		make_shared<Sphere>(float3(0, 25, -20), 1.f, nullptr),
		float3(1, .941, .914) * 200
	);
	// fill light
	scene->addAreaLight(
		make_shared<Sphere>(float3(10, 25, -20), 1.f, nullptr),
		float3(1, .941, .914) * 50
	);
	// back light
	scene->addAreaLight(
		make_shared<Sphere>(float3(0, 20, 10), 5.f, nullptr),
		float3(1, .941, .914) * 1
	);

	scene->camera.lookfrom = normalize(float3(1.24, 0.1, -4.84))*5;
	scene->camera.lookat = float3(0, 0, 0);
	scene->camera.vup = float3(0, 1, 0);
	scene->camera.aspect_ratio = 1;
	scene->camera.aperture = .1f;

	scene->lights.push_back(make_shared<InfiniteAreaLight>(float3(.765, .82, 1) * .05));
	//scene->lights.push_back(lightarea);
}

TheApp* CreateApp() { return new MyApp(); }

// -----------------------------------------------------------
// Initialize the application
// -----------------------------------------------------------
void MyApp::Init()
{
	scene = std::make_shared<Scene>();

	// anything that happens only once at application start goes here
	SimpleTestScene(scene.get());
	camera = make_shared<RotatingCamera>(scene->camera);

	integratorR = make_shared<PathTracer>();
	integratorL = integratorR; // make_shared<DbgIntegrator>();

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
	bool reportStats = true;
	// rotate camera if mouse moved since last tick
	if (!paused && (mouseDiff.x != 0 || mouseDiff.y != 0)) {
		float rotation_speed = 0.005f;
		camera->update(float2(-mouseDiff.y * rotation_speed, -mouseDiff.x * rotation_speed));
		accumulator->Clear();

		std::cerr << "camera.origin = " << camera->GetOrigin() << std::endl;
		reportStats = false;
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
				float3 clr = (uv.x < 0.5f) ? integratorL->Li(ray, *scene) : integratorR->Li(ray, *scene);
				if (HasNans(clr) || std::isinf(Luminance(clr))) {
					//std::cerr << "outlier\n";
					clr = float3(0.f);
				}
				accumulator->AddSample(x, y, clr);
			}
		}

		if (timer.stop() && reportStats) {
			auto stats = timer.getStatsAndReset();
			printf("Frame %d\tRender stats { best = %.2f, avg = %.2f, worst = %.2f }\t\t\r",
				accumulator->NumSamples(), stats.bestDuration, stats.avgDuration, stats.worstDuration);
		}

		if (numFramesToPause > 0 && accumulator->NumSamples() >= numFramesToPause) {
			numFramesToPause += 128;
			paused = true;
		}
	}


	accumulator->CopyToSurface(screen);

	if (paused && accumulator->IsInside(mousePos.x, mousePos.y)) {
		auto fc = accumulator->WindowToFilm(mousePos);
		Ray ray = camera->GetRay(fc.x, fc.y);

		if (leftMouseBtnDown) {
			integratorL->Li(ray, *scene);
		}

		SurfaceInteraction hit;
		if (scene->NearestIntersection(ray, hit)) {
			auto p = accumulator->FilmToWindow(camera->WorldToScreen(hit.p));
			auto n = accumulator->FilmToWindow(camera->WorldToScreen(hit.p + hit.shading.n*5));
			auto g = accumulator->FilmToWindow(camera->WorldToScreen(hit.p + hit.n * 5));
			auto dpdu = accumulator->FilmToWindow(camera->WorldToScreen(hit.p + hit.dpdu*5));
			auto dpdv = accumulator->FilmToWindow(camera->WorldToScreen(hit.p + hit.dpdv*5));
			screen->Line(p.x, p.y, n.x, n.y, 0x0000FF);
			screen->Line(p.x, p.y, g.x, g.y, 0xFFFF00);
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