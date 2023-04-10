#pragma once
#include <vector>

template <typename Predicate> 
int FindInterval(int size, const Predicate& pred) {
    int first = 0, len = size;
    while (len > 0) {
        int half = len >> 1, middle = first + half;
        // Bisect range based on value of _pred_ at _middle_
        if (pred(middle)) {
            first = middle + 1;
            len -= half + 1;
        } else
            len = half;
    }
    return clamp(first - 1, 0, size - 2);
}

struct Distribution1D {

    Distribution1D(const float* f, int n) : func(f, f + n), cdf(n + 1) {
        // Compute integral of step function at x_i
        cdf[0] = 0;
        for (int i = 1; i < n + 1; i++) cdf[i] = cdf[i - 1] + func[i - 1] / n;

        // Transform step function integral to CDF
        funcInt = cdf[n];
        if (funcInt == 0) {
            for (int i = 1; i < n + 1; i++) cdf[i] = float(i) / float(n);
        } else {
            for (int i = 1; i < n + 1; i++) cdf[i] /= funcInt;
        }
    }

    int Count() const { return (int)func.size(); }

    float SampleContinuous(float u, float* pdf, int* off = nullptr) const {
        // Find surrounding CDF segments and _offset_
        int offset = FindInterval((int)cdf.size(), [&](int index) { return cdf[index] <= u; });
        if (off) *off = offset;
        // Compute offset along CDF segment
        float du = u - cdf[offset];
        if ((cdf[offset + 1] - cdf[offset]) > 0) {
            du /= cdf[offset + 1] - cdf[offset];
        }

        // Compute PDF for sampled offset
        if (pdf) *pdf = (funcInt > 0) ? func[offset] / funcInt : 0;

        // Return x in range [0,1) corresponding to sample
        return (offset + du) / Count();
    }

    int SampleDiscrete(float u, float* pdf = nullptr, float* uRemapped = nullptr) const {
        // Find surrounding CDF segment and _offset_
        int offset = FindInterval((int)cdf.size(), [&](int index) {return cdf[index] <= u; });
        if (pdf) *pdf = (funcInt > 0) ? func[offset] / (funcInt * Count()) : 0;
        if (uRemapped)
            *uRemapped = (u - cdf[offset]) / (cdf[offset + 1] - cdf[offset]);
        return offset;
    }

    float DiscretePDF(int index) const {
        return func[index] / (funcInt * Count());
    }

    std::vector<float> func, cdf;
    float funcInt;
};