#pragma once

#include <vector>
#include <deque>
#include "argparser.h"
#include "readSplitcsv.h"
#include "gaussFilter.h"
#include <set>
#include "py_functions.h"
#include <cfloat>
#include <climits>
#include "ctpl_stl.h"
#include <unordered_set>
#include <stdint.h>
#include <math.h>
#include <assert.h>
#include <regex>
#include <array>


// below are some simple containers for holding data; 

// cigar string member
struct cig_element {
    char c;
    int count;
};

//All Information pertaining to an individual read
struct Read {
    int id;
    std::string name;
    std::string chr;
    char strand;
    int tint;
    std::vector<std::pair<std::array<int, 4>, std::vector<cig_element>>> intervals;
    std::string seq;
    int length;
    std::deque<uint8_t> data;
    std::vector<std::string> gaps;
};
struct Tint{
    int id;
    std::string chr;
    std::vector<std::array<int, 2>> intervals;
    int read_counts;
    std::vector<int> final_positions;
    std::vector<std::array<int, 2>> segs;
};

// struct holding data used by DP algorithm
// cov is for % of positions in genomic section ij
// OTher 3 members are the binary indicators
// in is for the number of reads with partial split 
// alignment coverage
struct mem{
    std::vector<double> cov;
    std::vector<bool> yay;
    std::vector<bool> nay;
    std::vector<bool> amb;
    int in;
};

// DP table used in DP algorithm
struct _DB {
    double D;
    std::array<int16_t, 3> B;
};


struct optimize_results {
    std::unordered_map<uint64_t, _DB> DB;
    double max_d;
    std::array<int16_t, 3> max_b;
    std::unordered_map<uint64_t, double> out_mem;
};

struct splicing_data {
    std::unordered_map<int, std::array<uint32_t, 2>> pos_to_Yy_idx;
    std::vector<std::vector<std::deque<bool>>> Yy_idx_to_r_indxs;
    std::vector<std::vector<uint32_t>> Yy_idx_to_pos;
    std::vector<std::vector<double>> Y_raw;
};


struct Ilp {
    int i;
    int l;
    double p;
    char c;
};


typedef  std::pair<Tint, std::vector<Read>> _tint_data_structure;


double round_mydouble(double var);

std::vector<double> smooth_threshold(double threshold);

_tint_data_structure read_split(std::string fp);

int forward_thread_cigar(const std::vector<cig_element> &cigar, int t_goal, int t_pos, int q_pos);

std::array<int, 2> get_interval_start(int start, const Read &read);

std::array<int, 2> get_interval_end(int end, const Read &read);

std::vector<Ilp> find_longest_poly_backwards(const std::string & seq, int s, int e, int step,
                                                 char _char='A', int match_score = 1, int mismatch_score = -2);

std::vector<Ilp> find_longest_poly(const std::string & seq, int s, int e, int step,
                                                 char _char='A', int match_score = 1, int mismatch_score = -2);

void get_unaligned_gaps_and_polyA(Read & read, const std::vector<std::array<int, 2>> & segs);

std::unordered_map<char, char>& getMap();

std::vector<std::vector<uint32_t>> get_cumulative_coverage(const std::vector<uint32_t> & candidate_y_indxs,
                         const std::vector<std::deque<bool>>& y_idx_to_r_idxs);

std::vector<std::vector<uint32_t>> get_cumulative_coverage(const std::deque<uint32_t> & candidate_y_indxs,
                         const std::vector<std::deque<bool>>& y_idx_to_r_idxs);

// The functions below are using simple bit shift operations to create keys
// for the DP algorithm. 
// I did this because Array Hashing is not provided by the STL library for 
// unordered_maps. I could have used boost or some other external library but 
// I wanted to avoid doing that to keep things simple.
// ___________________________________
// |XXXXXXXXXXXXXXXX|XXXXXXXXXXXXXXXX|
// <---------------><--------------->
// first entry = i   second entry = j  {both are 2 byte ints}
// second function follows a similar principle but with 3 int16_t's instead of 2
// As such the resulting key is a long instead of a 4byte integer

uint32_t _2array_to_idx(const std::array<int16_t, 2>&& arr);
uint64_t _3array_to_idx(const std::array<int16_t, 3>&& arr);

uint32_t _2array_to_idx(const std::array<int16_t, 2>& arr);
uint64_t _3array_to_idx(const std::array<int16_t, 3>& arr);

std::vector<int> run_optimize(const std::deque<uint32_t> & candidate_y_idxs, const std::vector<int> & fixed_c_idxs,
                  const std::vector<std::vector<uint32_t>> & coverage,
                  const std::vector<double> & smoothed_thrshold, double threshold_rate,
                  int min_read_support_outside);

optimize_results optimize(const std::deque<uint32_t> & candidate_y_idxs,
              const std::vector<std::vector<uint32_t>> & C,
              int start, int end, const std::vector<double>& smoothed_threshold, double threshold_rate,
              int read_support);

double get_high_threshold(uint32_t seg_len, const std::vector<double> &smoothed_threshold, double threshold_rate);

void read_sequence(_tint_data_structure &tint, const std::string & fp);

std::deque<uint32_t> candidates_from_peaks(std::vector<double> & y, uint32_t f_y_idx, uint32_t l_y_idx);

void break_large_problems(const std::deque<uint32_t> & candidate_y_idxs,
                                        std::vector<int> & fixed_c_idxs, const std::vector<double> & y,
                                        uint32_t max_problem_size, int & new_problems_total_count, int window = 5);


void refine_segmentation(const std::vector<double> &y_raw, std::vector<uint32_t> &y_idxs,
                                    double sigma, int skip = 20, int min_interval_splice = 20);
splicing_data process_splicing_data(const _tint_data_structure & tint);

void run_segment(int id, const segArgs & seg);

// To keep things from running too slow I impletemented the sets from the original code as
// vectors, running the function below to update the set when necessary
void fast_set(std::vector<int>& v);

void fast_set(std::vector<std::string>& v);

int segment(_tint_data_structure &tint, double sigma, const std::vector<double> &smoothed_threshold,
             double threshold_rate, double variance_factor, uint32_t max_problem_size, int min_read_support_outside);
