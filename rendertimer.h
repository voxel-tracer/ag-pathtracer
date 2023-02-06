#pragma once

struct FrameStats {
	float bestDuration;
	float avgDuration;
	float worstDuration;
};

class RenderTimer {
public:
	void start() {
		timer.reset();
	}

	// returns true if it is time to print frame stats
	bool stop() {
		float frameDuration = timer.elapsed() * 1000;
		frames++;
		duration += frameDuration;
		if (frameDuration < bestDuration)
			bestDuration = frameDuration;
		if (frameDuration > worstDuration)
			worstDuration = frameDuration;

		return duration >= sampleDuration;
	}

	FrameStats getStatsAndReset() {
		FrameStats stats{
			bestDuration,
			duration / frames,
			worstDuration
		};

		frames = 0;
		duration -= sampleDuration;
		bestDuration = FLT_MAX;
		worstDuration = 0;

		return stats;
	}

private:
	Timer timer;
	float sampleDuration = 1000.f; // time interval used to compute best/worst durations

	int frames = 0;
	float duration = 0;
	float bestDuration = FLT_MAX;
	float worstDuration = 0;
};
