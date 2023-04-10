#pragma once

class Texture {
public:
	virtual float3 value(float u, float v) const = 0;
};

class SolidColor : public Texture {
public:
	SolidColor(float a) : SolidColor(a, a, a) {}
	SolidColor(float r, float g, float b) : SolidColor(float3(r, g, b)) {}
	SolidColor(float3 c) : color(c) {}

	virtual float3 value(float u, float v) const override {
		return color;
	}

	static shared_ptr<SolidColor> make(float3 v) {
		return make_shared<SolidColor>(v);
	}

private:
	float3 color;
};

class CheckerTexture : public Texture {
public:
	CheckerTexture(shared_ptr<Texture> c1, shared_ptr<Texture> c2) : color1(c1), color2(c2) {}

	virtual float3 value(float u, float v) const override {
		if (((int)(u * 8) + (int)(v * 8)) % 2)
			return color1->value(u, v);
		return color2->value(u, v);
	}

private:
	shared_ptr<Texture> color1;
	shared_ptr<Texture> color2;
};

class HDRTexture : public Texture {
public:
	HDRTexture(const std::string& filename) {
		int n;
		float* data = stbi_loadf(filename.c_str(), &width, &height, &n, 0);
		if (data) {
			pixels = (float3*)MALLOC64(width * height * sizeof(float3));
			const int s = width * height;
			for (int i = 0; i < s; i++) 
				pixels[i] = float3(data[i * n + 0], data[i * n + 1], data[i * n + 2]);
		}
		stbi_image_free(data);
	}

	~HDRTexture() {
		FREE64(pixels);
	}

	float3 value(float u, float v) const {
		int s = (int)std::floor(u * width - .5f);
		int t = (int)std::floor(v * height - .5f);

		// just clamp the values if they are outside [0,1] range
		int x = Mod(s, width);
		int y = Mod(t, height);
		return GetPixel(x, y);
	}

	int Width() const { return width; }
	int Height() const { return height; }

	float3 GetPixel(int x, int y) const {
		return pixels[y * width + x];
	}

	static int Mod(int a, int b) {
		int result = a - (a / b) * b;
		return (int)((result < 0) ? result + b : result);
	}

private:
	int width, height;
	float3* pixels;
};
