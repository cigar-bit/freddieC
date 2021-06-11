#pragma once

#include <cmath>
#include <vector>
#include <assert.h>
#include <iostream>
#include <numeric>
#include <algorithm>


// This header as well as the py_functions header provide some simple c++ implementations of 
// the gaussian_filter1d and find_peaks functions from scipy. 

inline std::vector<double> _gaussian_kernel1d(double sigma, int radius) {
    double sigma_sqr = sigma*sigma;
    std::vector<double> phi_x(radius*2 + 1);
    std::vector<int> x(radius*2 + 1);
    std::iota(x.begin(), x.end(), -radius);

    for (uint32_t i = 0; i < phi_x.size(); i++) {
        phi_x[i] = std::exp(-0.5f / sigma_sqr * x[i] * x[i]);
    }
    double sum = std::accumulate(phi_x.begin(), phi_x.end(), 0.0f);
    for (auto & p: phi_x) {
        p /= sum;
    }
    return phi_x;
}


inline std::vector<double> dirty_gaussian_filter1d(std::vector<double> & input, double sigma,
                                              double truncate = 4.0f, bool constant = false) {
    int lw = truncate * sigma + 0.5f;
    std::vector<double> weights = _gaussian_kernel1d(sigma, lw);
    std::reverse(weights.begin(), weights.end());
    std::vector<double> output(input.size());
    int idx_shift = weights.size() / 2;
    int count;
    if (lw < input.size() && !constant){

        for (int i = 0; i < (int)input.size(); i++) {
            count = i - idx_shift;
            for (int j = 0; j < (int)weights.size(); j++) {
                if (count < 0) output[i] += input[-count - 1]*weights[j];
                else if (count >= (int)input.size()) output[i] += 
                        input[2*input.size() - count - 1]*weights[j];
                else output[i] += input[count]*weights[j];
                count++;
            }
        }
    }
    else {
        for (int i = 0; i < (int)input.size(); i++) {
            count = i - idx_shift;
            for (int j = 0; j < (int)weights.size(); j++) {
                if (count < 0) output[i] += input.front()*weights[j];
                else if (count >= (int)output.size()) output[i] += input.back()*weights[j];
                else output[i] += input[count]*weights[j];
                count++;
            }
        }
    }
    return output;
}
