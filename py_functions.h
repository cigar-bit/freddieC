#pragma once
#include <deque>
#include <stdint.h>
#include <vector>
#include <limits>
#include <algorithm>
#include <tgmath.h>

// This header as well as the gaussFilter header provide some simple c++ implementations of 
// the gaussian_filter1d and find_peaks functions from scipy. 

inline bool AreSame(double a, double b)
{
    return std::fabs(a - b) < std::numeric_limits<double>::epsilon();
}

template <class Iter>
std::deque<uint32_t> find_peaks(Iter a, Iter b, Iter c, const std::vector<double> &v, int distance = 0.0f) {
    std::deque<uint32_t> temp; 
    a++;
    int rshift = 1; 
    while(a != b - 1) {
        rshift = 1;
        while (*a == *(a + rshift) && a + rshift != b) rshift++;
        if (*a > *(a - 1) && *a > *(a + rshift)) temp.push_back(a - c);
        a++;
    }
    double max = -std::numeric_limits<double>::max(); uint32_t max_idx;
    for (auto i : temp) {
        if (max < v.at(i)) {
            max = v[i]; max_idx = i;
        }
    }

    if (distance != 0) {
        for (auto it = temp.begin(); it != temp.end();) {
            bool del_peak;
            if (it == temp.begin()) del_peak = ( (*(it + 1) - *it < distance) && 
                                            (v.at(*it) < v.at(*(it + 1))));
            else if (it == temp.end() - 1) del_peak = ((*it - *(it - 1) < distance) && 
                                                (v.at(*it) < v.at(*(it - 1))));
            else del_peak = ((*it - *(it - 1) < distance) && (v.at(*it) < v.at(*(it - 1)))) || 
                        ((*(it + 1) - *it < distance) && (v.at(*it) < v.at(*(it + 1))));
            if (del_peak) it = temp.erase(it);
            else ++it;
        }
    }

    return temp;
}
