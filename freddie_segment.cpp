#include "freddie_segment.h"

static const std::string cinterval_re = "[0-9]+-[0-9]+";
static const std::regex celement_re("-");
static const std::regex comma_prog(",");
static const std::regex colon_prog(":");
static const std::regex relement_re("-");
static const std::regex tab_re(R"(\t)");
static const std::string cigar_op_re = "[0-9]+[MIDNSHPX=]";
static const std::array<char, 2> char_holder = {'A', 'T'};
static const std::string cigar_re = [](){
    char c[40];
    std::string s;
    sprintf(c, "(?:%s)+", cigar_op_re.c_str());
    s = c;
    return s;
}();

static const std::string rinterval_re = [](){
    char c[200];
    std::string s = "[0-9]+-[0-9]+";
    sprintf(c, "(?:%s:%s:%s)+", s.c_str(), s.c_str(), cigar_re.c_str());
    s = c;
    return s;
}();

static const std::string tint_prog_sub = [](){
    char c[200];
    std::string s;
    sprintf(c, "%s(?:,(?:%s))*", cinterval_re.c_str(), cinterval_re.c_str());
    s = c;
    return s;
}();

static const std::string read_prog_sub = [](){
    char c[400];
    std::string s;
    sprintf(c, "%s(:?\\t(?:%s))*", rinterval_re.c_str(), rinterval_re.c_str());
    s = c;
    return s;
}();

static const std::array<std::string, 4> tint_prog = {std::string(R"([0-9A-Za-z!#$%&+./:;?@^_|~-][0-9A-Za-z!#$%&*+./:;=?@^_|~-]*)"),
                                                    std::string(R"([0-9]+)"),
                                                    std::string(tint_prog_sub),
                                                    std::string(R"([0-9]+)")
                                                   };
static const std::regex tint_prog_re = [](){
    char c[800];
    std::string s;
    sprintf(c, "#(%s)\\t(%s)\\t(%s)\\t(%s)", tint_prog[0].c_str(), tint_prog[1].c_str(),
                                tint_prog[2].c_str(), tint_prog[3].c_str());
    s = c;
    return std::regex(s);
}();

// tint prog index layout:
// 0 chr_re
// 1 cid_re
// 2 intervals_re
// 3 read_count_re

static const std::array<std::string, 6> read_prog = {std::string(R"([0-9]+)"),
                                                    std::string(R"([!-?A-~]{1,254})"),
                                                    std::string(R"([0-9A-Za-z!#$%&+./:;?@^_|~-][0-9A-Za-z!#$%&*+./:;=?@^_|~-]*)"),
                                                    std::string(R"([+-])"),
                                                    std::string(R"([0-9]+)"),
                                                    std::string(read_prog_sub)
                                                   };

static const std::regex read_prog_re = [](){
    char c[800];
    std::string s;
    sprintf(c, R"((%s)\t(%s)\t(%s)\t(%s)\t(%s)\t(%s))", read_prog[0].c_str(), read_prog[1].c_str(),
                                                        read_prog[2].c_str(), read_prog[3].c_str(),
                                                        read_prog[4].c_str(), read_prog[5].c_str());
    s = c;
    return std::regex(s);
}();

// read prog index layout:
// 0 rid_re
// 1 name_re
// 2 chr_re
// 3 strand_re
// 4 cid_re
// 5 intervals_re

static const std::regex cinterval_prog(R"(-)");
static const std::regex rinterval_prog(R"(-)");
static const std::regex cigar_prog(R"([0-9]+[MIDNSHPX=])");


std::unordered_map<char, char>& getMap() {
    static std::unordered_map<char, char> rev_comp{ {'A', 'T'}, {'T', 'A'}, {'C', 'G'}, {'G', 'C'} };
    return rev_comp;
}




double round_mydouble(double var, int x) {
    x = std::pow(10, x) + 0.5f;
    double value = (int)(var * x + 0.5f);
    return (double)value / x;
}

std::vector<double> smooth_threshold(double threshold) {
    std::vector<double> smooth; smooth.reserve(1000);
    uint32_t x;
    double y;
    while (true) {
        x = smooth.size();
        y = threshold/(1 + ((threshold-0.5f)/0.5f)*exp(-0.05*x));
        if (x > 5 && x*(threshold-y) < 0.5) break;
        smooth.emplace_back(round_mydouble(y, 2));
        assert(y < 1000);
    }
    smooth.shrink_to_fit();
    return smooth;
}

_tint_data_structure read_split(std::string fp) {
    tsv_table_delineated split_tsv = read_tsv_delineated(fp);
    _tint_data_structure tints;
    for (auto & line : split_tsv) {
        if (line.at(0) == '#') {
            std::smatch res;
            std::regex_search(line, res, tint_prog_re);
            std::deque<std::string> matches;
            const std::string & s = res[3].str();
            std::copy(std::sregex_token_iterator(s.begin(),
                      s.end(), comma_prog, -1),
                    std::sregex_token_iterator(),
                    std::back_insert_iterator(matches));

            std::array<std::string, 2> str_matcher;
            std::vector<std::array<int, 2>> intervals(matches.size());
            for (uint32_t i = 0; i < matches.size(); i++) {
                std::copy(std::sregex_token_iterator(matches[i].begin(),
                          matches[i].end(), celement_re, -1),
                        std::sregex_token_iterator(),
                        std::begin(str_matcher));

                intervals[i] = {std::stoi(str_matcher[0]), std::stoi(str_matcher[1])};
            }
            Tint tint = {std::stoi(res[2].str()), res[1].str(), intervals,
                         std::stoi(res[4].str()), {}};
            tints.first = std::move(tint);
        }
        else {
            std::smatch res;

            std::regex_search(line, res, read_prog_re);

            const std::string & full_s = res[6].str();
            std::deque<std::string> all_matches;
            std::copy(std::sregex_token_iterator(full_s.begin(),
                      full_s.end(), tab_re, -1),
                    std::sregex_token_iterator(),
                    std::back_insert_iterator(all_matches));
            std::vector<std::pair<std::array<int, 4>, std::vector<cig_element>>> intervals(all_matches.size());
            int count = 0;
            for (auto & am : all_matches) {
                std::deque<std::string> matches;
                const std::string & s = am;
                std::copy(std::sregex_token_iterator(s.begin(),
                          s.end(), colon_prog, -1),
                        std::sregex_token_iterator(),
                        std::back_insert_iterator(matches));


                for (int i = 0; i < 2; i++) {
                    std::array<std::string, 4> c_matches;
                    std::string & c_s = matches[i];
                    std::copy(std::sregex_token_iterator(c_s.begin(),
                              c_s.end(), celement_re, -1),
                            std::sregex_token_iterator(),
                            std::begin(c_matches));
                    intervals[count].first[i*2] = std::stoi(c_matches[0]);
                    intervals[count].first[i*2 + 1] = std::stoi(c_matches[1]);
                }
                std::deque<std::string> c_matches;
                std::string & c_s = matches[2];
                std::copy(std::sregex_token_iterator(c_s.begin(),
                          c_s.end(), cigar_prog, 0),
                        std::sregex_token_iterator(),
                        std::back_insert_iterator(c_matches));
                intervals[count].second.reserve(c_matches.size());
                for (auto& m : c_matches) {
                    char c = m.back();
                    m.pop_back();
                    intervals[count].second.push_back({c, std::stoi(m)});
                }
                count++;
            }


            Read read = {std::stoi(res[1].str()), res[2].str(), res[3].str(),
                         res[4].str()[0], std::stoi(res[5].str()), intervals, "", 0};
            tints.second.emplace_back(read);
        }
    }
    return tints;
}


void read_sequence(_tint_data_structure & tint, const std::string & fp)
{
    tsv_table table = read_tsv(fp);
    std::unordered_map<int, std::string> rid_to_seq;
    for (auto& line : table) {
        rid_to_seq[std::stoi(line[0])] = std::move(line[3]);
    }
    for (auto & read : tint.second) {
        read.seq = rid_to_seq.at(read.id);
        read.length = rid_to_seq.at(read.id).size();
    }
}



void run_segment(int id, const segArgs & seg) {

    _tint_data_structure tint = read_split(seg.split_dir + "/" + seg.contig + "/split_" + seg.contig +
                                           "_" + std::to_string(seg.tint_id) + ".tsv");

    read_sequence(tint, seg.split_dir + "/" + seg.contig + "/reads_" + seg.contig +
                  "_" + std::to_string(seg.tint_id) + ".tsv");
    segment(tint, seg.sigma, seg.smooth_threshold, seg.threshold_rate, seg.variance_factor,
            seg.max_problem_size, seg.min_read_support_outside);
    std::string outdir_fp = seg.out_dir + seg.contig + "/segment_" + seg.contig + "_" +
            std::to_string(tint.first.id) + ".tsv";
    tsv_table record(tint.second.size() + 1, std::vector<std::string>(7));
    record[0].resize(3);
    record[0][0] = "#" + tint.first.chr;
    record[0][1] = std::to_string(tint.first.id);
    record[0][2] = join_container(tint.first.final_positions, ",");
    int idx = 1;
    for (auto & read : tint.second) {
        record[idx][0] = std::to_string(read.id);
        record[idx][1] = read.name;
        record[idx][2] = read.chr;
        record[idx][3].push_back(read.strand);
        record[idx][4] = std::to_string(read.tint);
        record[idx][5] = join_container(read.data, ""); 
        if (!read.gaps.empty()) record[idx][6] = join_container(read.gaps, ",") + ",";
        idx++;
    }
    write_tsv_table(record, outdir_fp);
}

void fast_set(std::vector<int>& v) {
    std::unordered_set<int> s;
    for (int i : v) s.insert(i);
    v.assign(s.begin(), s.end());
    std::sort(v.begin(), v.end());
}

void fast_set(std::vector<std::string>& v) {
    std::unordered_set<std::string> s;
    for (auto & e : v) s.insert(e);
    v.assign(s.begin(), s.end());
    std::sort(v.begin(), v.end());
}


int segment(_tint_data_structure & tint, double sigma, const std::vector<double> & smoothed_threshold,
             double threshold_rate, double variance_factor, uint32_t max_problem_size,
             int min_read_support_outside)
{
    splicing_data splice = process_splicing_data(tint);
    std::vector<std::vector<double>> Y(splice.Y_raw.size()); 
    for (uint32_t i = 0; i < Y.size(); i++) {
        Y[i] = dirty_gaussian_filter1d(splice.Y_raw[i], sigma, 4.0f);
    }
    std::deque<double> Y_none_zero_vals; 
    for (auto & y : Y) {
        for (auto & v : y) {
            if (v > 0) Y_none_zero_vals.push_back(v);
        }
    }
    double mean = std::accumulate(Y_none_zero_vals.begin(),
                                  Y_none_zero_vals.end(), 0.0f) / Y_none_zero_vals.size();

    std::vector<double> diff(Y_none_zero_vals.size());
    std::transform(Y_none_zero_vals.begin(), Y_none_zero_vals.end(),
                   diff.begin(), [mean](double x) { return x - mean; });
    double sq_sum = std::inner_product(diff.begin(), diff.end(), diff.begin(), 0.0);
    double stdev = std::sqrt(sq_sum / Y_none_zero_vals.size());


    double variance_threshold = mean + variance_factor * stdev;
    int Y_idx = 0;
    for (auto & y : Y) {
        std::vector<int> fixed_c_idxs; fixed_c_idxs.reserve(1000);

        uint32_t f_y_idx = 0;
        uint32_t l_y_idx = y.size() - 1;
        fixed_c_idxs.push_back(0);

        std::deque<uint32_t> current_candidates = candidates_from_peaks(y, f_y_idx, l_y_idx);

        current_candidates.push_front(f_y_idx);
        current_candidates.push_back(l_y_idx);
        std::sort(current_candidates.begin(), current_candidates.end());
        fixed_c_idxs.push_back(current_candidates.size() - 1);
        int c_idx = 0;
        for (auto & y_idx : current_candidates) {
            if (y[y_idx] > variance_threshold) fixed_c_idxs.push_back(c_idx);
            c_idx++;
        }

        int new_problems_total_count;
        fast_set(fixed_c_idxs);
        break_large_problems(current_candidates, fixed_c_idxs, y, max_problem_size,
                             new_problems_total_count);


        std::vector<std::vector<uint32_t>> cumulative_coverage =
                get_cumulative_coverage(current_candidates, splice.Yy_idx_to_r_indxs[Y_idx]);
        std::sort(fixed_c_idxs.begin(), fixed_c_idxs.end());
        std::vector<int> final_c_idxs = run_optimize(current_candidates, fixed_c_idxs, cumulative_coverage,
                                                  smoothed_threshold, threshold_rate,
                                                  min_read_support_outside);
        std::vector<uint32_t> final_y_idxs; final_y_idxs.reserve(final_c_idxs.size());
        for (auto & c_idx : final_c_idxs) {
            final_y_idxs.push_back(current_candidates[c_idx]);
        }
        refine_segmentation(splice.Y_raw[Y_idx], final_y_idxs, sigma);
        std::sort(final_y_idxs.begin(), final_y_idxs.end());
        int y_idx_counter = tint.first.final_positions.size();
        tint.first.final_positions.resize(tint.first.final_positions.size() + final_y_idxs.size());
        for (auto & y_idx : final_y_idxs) {
            tint.first.final_positions[y_idx_counter] = splice.Yy_idx_to_pos[Y_idx][y_idx];
            y_idx_counter++;
        }
        cumulative_coverage = get_cumulative_coverage(final_y_idxs, splice.Yy_idx_to_r_indxs[Y_idx]);
        int seg_idx = 0;
        int seg_len;
        for (std::vector<uint32_t>::iterator it = final_y_idxs.begin(); it != --final_y_idxs.end();) {
            int s_yidx = *it;
            int e_yidx = *(++it);
            seg_len = splice.Yy_idx_to_pos[Y_idx][e_yidx] - splice.Yy_idx_to_pos[Y_idx][s_yidx]+1;
            double h = get_high_threshold(seg_len, smoothed_threshold, threshold_rate);
            double l = 1-h;
            int r_idx = 0;
            for (auto & read : tint.second) {
                double cov_ratio = (double)(cumulative_coverage[seg_idx+1][r_idx] -
                        cumulative_coverage[seg_idx][r_idx])/seg_len;
                if (cov_ratio > h) read.data.push_back(1);
                else if (cov_ratio < l) read.data.push_back(0);
                else read.data.push_back(2);
                r_idx++;
            }
            seg_idx++;
        }
        for (auto & read : tint.second) {
            read.data.push_back(0);
        }
        Y_idx++;
    }

    tint.first.segs.reserve(tint.first.final_positions.size() - 1);
    for (std::vector<int>::iterator it = tint.first.final_positions.begin();
         it != tint.first.final_positions.end() - 1;) {
        int s = *it;
        int e = *(++it);
        tint.first.segs.push_back({s, e});
    }
    for (auto & read : tint.second) {
        read.data.pop_back();
        get_unaligned_gaps_and_polyA(read, tint.first.segs);
    }
    return tint.first.id;
}

splicing_data process_splicing_data(const _tint_data_structure & tint)
{
    splicing_data splice = {
        std::unordered_map<int, std::array<uint32_t, 2>>(),
        std::vector<std::vector<std::deque<bool>>>(tint.first.intervals.size()),
        std::vector<std::vector<uint32_t>>(tint.first.intervals.size()),
        std::vector<std::vector<double>>(tint.first.intervals.size())
    };
    uint32_t count = 0; 
    for (auto & t: tint.first.intervals) {
        std::vector<uint32_t> y_idx_to_pos(t[1] - t[0] + 1);
        for (uint32_t p = t[0]; p < t[1] + 1; p++) {
            splice.pos_to_Yy_idx.insert(std::pair<int, std::array<uint32_t, 2>>(p, {count, p - t[0]})); 
            y_idx_to_pos[p - t[0]] = p;
        }
        splice.Yy_idx_to_r_indxs[count] = std::vector<std::deque<bool>>
                                        (y_idx_to_pos.size(),
                                            std::deque<bool>(tint.second.size()));
        splice.Y_raw[count] = std::vector<double>(y_idx_to_pos.size());
        splice.Yy_idx_to_pos[count] = y_idx_to_pos;
        count++;
    }
    int r_idx = 0;
    for (auto & read : tint.second) {
        for (auto & in : read.intervals) {
            unsigned int & Y_idx_s = splice.pos_to_Yy_idx.at(in.first[0])[0];
            unsigned int & y_idx_s = splice.pos_to_Yy_idx.at(in.first[0])[1];
            unsigned int & y_idx_e = splice.pos_to_Yy_idx.at(in.first[1])[1];
            splice.Y_raw.at(Y_idx_s).at(y_idx_s) += 1.0f;
            splice.Y_raw.at(Y_idx_s).at(y_idx_e) += 1.0f;
            auto & bool_vec = splice.Yy_idx_to_r_indxs.at(Y_idx_s);
            for (uint32_t i = y_idx_s; i < y_idx_e; i++) bool_vec[i][r_idx] = true;
        }
        r_idx++;
    }
    return splice;
}

std::deque<uint32_t> candidates_from_peaks(std::vector<double> &y, uint32_t f_y_idx, uint32_t l_y_idx)
{
    std::deque<uint32_t> c = find_peaks(y.begin() + f_y_idx, y.begin() + l_y_idx + 1, y.begin(), y);
    for (auto & e : c) e += f_y_idx;
    return c;
}



void break_large_problems(const std::deque<uint32_t> &candidate_y_idxs,
                          std::vector<int> &fixed_c_idxs, const std::vector<double> &y,
                          uint32_t max_problem_size, int &new_problems_total_count, int window)
{
    std::vector<std::array<int, 2>> fixed_c_idxs_pairs(fixed_c_idxs.size() - 1);

    int count = 0;
    for (std::vector<int>::iterator it = fixed_c_idxs.begin(); it != --fixed_c_idxs.end(); ++it) {
        fixed_c_idxs_pairs[count] = { *it, *(++it)--};
        count++;
    }
    new_problems_total_count = 0;
    uint32_t problem_size;
    for (const auto & arr : fixed_c_idxs_pairs) {
        const uint32_t & c_idx_s = arr[0];
        const uint32_t & c_idx_e = arr[1];
        problem_size = c_idx_e - c_idx_s + 1;
        if (problem_size <= max_problem_size) continue;
        uint32_t new_problems_count = (problem_size + max_problem_size - 1) / max_problem_size;
        uint32_t new_problems_size = problem_size / new_problems_count;
        new_problems_total_count += new_problems_count - 1;
        double max_c_idx_y_v; 
        int max_c_idx; 
        int mid_anchor;
        for (uint32_t i = 1; i < new_problems_count; i++) {
            mid_anchor = c_idx_s + i * new_problems_size;
            max_c_idx_y_v = -std::numeric_limits<double>::max();
            max_c_idx = -1;
            for (uint32_t c_idx = mid_anchor - window; c_idx < mid_anchor + window; c_idx++) {
                if (y[candidate_y_idxs[c_idx]] > max_c_idx_y_v) {
                    max_c_idx_y_v = y[candidate_y_idxs[c_idx]];
                    max_c_idx = c_idx;
                }
            }
            if (max_c_idx != -1) fixed_c_idxs.push_back(max_c_idx);
        }
    }
}

std::vector<std::vector<uint32_t>> get_cumulative_coverage(const std::vector<uint32_t> &candidate_y_indxs,
                                                       const std::vector<std::deque<bool>>& y_idx_to_r_idxs)
{
    std::vector<std::vector<uint32_t>> C(candidate_y_indxs.size() + 1,
                                         std::vector<uint32_t>(y_idx_to_r_idxs[0].size()));
    uint32_t C_idx = 1;

    std::vector<std::array<uint32_t, 2>> zipped_pairs(candidate_y_indxs.size() - 1);
    int count = 0;
    for (std::vector<uint32_t>::const_iterator it = candidate_y_indxs.begin();
         it != --candidate_y_indxs.end();)
    {
        zipped_pairs[count] = { *it, *(++it)};
        count++;
    }
    for (const auto & p : zipped_pairs) {
        const uint32_t & cur_y_idx = p[0];
        const uint32_t & nxt_y_idx = p[1];
        for (uint32_t i = 0; i < C[C_idx].size(); i++) {
            for (uint32_t j = cur_y_idx; j < nxt_y_idx; j++) {
                C[C_idx][i] += y_idx_to_r_idxs[j][i];
            }
        }
        C_idx++;
    }
    for (C_idx = 1; C_idx < C.size();C_idx++)  {
        for (uint32_t j = 0; j < C[0].size(); j++) {
            C[C_idx][j] += C[C_idx - 1][j];
        }
    }
    return C;
}

std::vector<std::vector<uint32_t>> get_cumulative_coverage(const std::deque<uint32_t> &candidate_y_indxs,
                                                       const std::vector<std::deque<bool>>& y_idx_to_r_idxs)
{
    std::vector<std::vector<uint32_t>> C(candidate_y_indxs.size() + 1,
                                         std::vector<uint32_t>(y_idx_to_r_idxs[0].size()));
    uint32_t C_idx = 1;

    std::vector<std::array<uint32_t, 2>> zipped_pairs(candidate_y_indxs.size() - 1);
    int count = 0;
    for (std::deque<uint32_t>::const_iterator it = candidate_y_indxs.begin();
         it != --candidate_y_indxs.end();)
    {
        zipped_pairs[count] = { *it, *(++it)};
        count++;
    }
    for (const auto & p : zipped_pairs) {
        const uint32_t & cur_y_idx = p[0];
        const uint32_t & nxt_y_idx = p[1];
        for (uint32_t i = 0; i < C[C_idx].size(); i++) {
            for (uint32_t j = cur_y_idx; j < nxt_y_idx; j++) {
                C[C_idx][i] += y_idx_to_r_idxs[j][i];
            }
        }
        C_idx++;
    }
    for (C_idx = 1; C_idx < C.size();C_idx++)  {
        for (uint32_t j = 0; j < C[0].size(); j++) {
            C[C_idx][j] += C[C_idx - 1][j];
        }
    }
    return C;
}


std::vector<int> run_optimize(const std::deque<uint32_t> &candidate_y_idxs, const std::vector<int> &fixed_c_idxs,
                  const std::vector<std::vector<uint32_t> > &coverage,
                  const std::vector<double> &smoothed_thrshold, double threshold_rate,
                  int min_read_support_outside)
{
    std::vector<int> final_c_idxs = fixed_c_idxs;
    uint32_t start, end;
    std::array<int16_t, 3> end_max_b = {-1, -1, -1};
    for (std::vector<int>::const_iterator it = fixed_c_idxs.begin(); it != --fixed_c_idxs.end();) {
        start = *it;
        end = *(++it);
        optimize_results res = optimize(candidate_y_idxs, coverage, start, end, smoothed_thrshold,
                                        threshold_rate, min_read_support_outside);
        while (res.max_b != end_max_b) {
            for (auto & e : res.max_b) final_c_idxs.push_back(e);
            res.max_b = res.DB.at(_3array_to_idx(res.max_b)).B;
        }
    }
    fast_set(final_c_idxs);
    return final_c_idxs;
}

uint32_t _2array_to_idx(const std::array<int16_t, 2>&& arr) {
    return ((int32_t)arr[0] << 16) | arr[1];
}


uint64_t _3array_to_idx(const std::array<int16_t, 3>&& arr) {
    return (((int64_t)arr[0] << 32) | ((int64_t)arr[1] << 16) | arr[2]);
}

uint32_t _2array_to_idx(const std::array<int16_t, 2>& arr) {
    return ((int32_t)arr[0] << 16) | arr[1];
}


uint64_t _3array_to_idx(const std::array<int16_t, 3>& arr) {
    return (((int64_t)arr[0] << 32) | ((int64_t)arr[1] << 16) | arr[2]);
}

optimize_results optimize(const std::deque<uint32_t> &candidate_y_idxs,
                            const std::vector<std::vector<uint32_t> > &C, 
                            int start, int end, const std::vector<double>& smoothed_threshold, 
                            double threshold_rate, int read_support)
{

    std::unordered_map<uint32_t, mem> mem_map;
    int seg_len;
    for (int16_t i = start; i < end; i++) {
        for (int16_t j = i; j < end + 1; j++) { 
            seg_len = candidate_y_idxs[j] - candidate_y_idxs[i] + 1;
            std::vector<double> temp_cov(C[0].size());
            for (int k = 0; k < (int)temp_cov.size(); k++) temp_cov[k] =
                    (double)(C[j][k] - C[i][k]) / seg_len;
            double h = get_high_threshold(seg_len, smoothed_threshold, threshold_rate);
            double l = 1 - h;
            std::vector<bool> temp_yay(C[0].size());
            for (int k = 0; k < temp_cov.size(); k++) temp_yay[k] =
                    temp_cov[k] > h;
            std::vector<bool> temp_nay(C[0].size());
            for (int k = 0; k < temp_cov.size(); k++) temp_nay[k] =
                    temp_cov[k] < l;
            std::vector<bool> temp_amb(C[0].size());
            for (int k = 0; k < temp_cov.size(); k++) temp_amb[k] =
                    !(temp_yay[k] || temp_nay[k]);

            mem_map[_2array_to_idx({i, j})] = {std::move(temp_cov), std::move(temp_yay), std::move(temp_nay),
                                std::move(temp_amb), INT_MAX};
        }
 
    }
    int mem_size = mem_map.size();
    auto inside = [&](const std::array<int16_t, 2> && ij){
        uint32_t ij_idx = _2array_to_idx(ij);
        int in_mem_size = mem_map.size();
        if (mem_map.at(ij_idx).in == INT_MAX) {
            if (ij[0] == ij[1])mem_map.at(ij_idx).in = 0;
            else {
                mem_map.at(ij_idx).in = 0;
                for (const auto & e : mem_map.at(ij_idx).amb) mem_map.at(ij_idx).in += e;
                mem_map.at(ij_idx).in *= -1;
            }
        }
        return mem_map.at(ij_idx).in;
    };
    std::unordered_map<uint64_t, double> out_mem;
    uint32_t idx_buffer[2];
    auto outside = [&](const std::array<int16_t, 3>& ijk) {
        int out_mem_size = out_mem.size();
        uint64_t ijk_idx = _3array_to_idx(ijk);
        if (out_mem.count(ijk_idx) == 0) {
            if (ijk[0] == ijk[1] || ijk[1] == ijk[2]) out_mem[ijk_idx] = 0;
            else {
                out_mem[ijk_idx] = 0;
                idx_buffer[0] = _2array_to_idx({ijk[0], ijk[1]});
                idx_buffer[1] = _2array_to_idx({ijk[1], ijk[2]});
                for (int32_t i = 0; i < mem_map.begin()->second.yay.size(); i++) {
                    out_mem[ijk_idx] += (mem_map.at(idx_buffer[0]).yay[i] && 
                    mem_map.at(idx_buffer[1]).nay[i]) ||
                    (mem_map.at(idx_buffer[1]).yay[i] && 
                    mem_map.at(idx_buffer[0]).nay[i]);
                }
            }
            if (out_mem.at(ijk_idx) < read_support) out_mem.at(ijk_idx) = -std::numeric_limits<double>::max();
        }
        return out_mem.at(ijk_idx);
    };
    std::unordered_map<uint64_t, _DB> DB;
    std::function<double(const std::array<int16_t, 3>&)> dp = [&](const std::array<int16_t, 3>& ijk) {
        int32_t watchsize1 = DB.size();
        uint64_t ijk_idx = _3array_to_idx(ijk);
        if (DB.count(ijk_idx) > 0) {
            return DB.at(ijk_idx).D;
        }
        std::array<int16_t, 3> dp_max_b = {-1, -1, -1};
        double dp_max_d = -std::numeric_limits<double>::max();
        if (candidate_y_idxs[ijk[1]]-candidate_y_idxs[ijk[0]] < 5 or
                candidate_y_idxs[ijk[2]]-candidate_y_idxs[ijk[1]] < 5) {
            DB[ijk_idx] = {dp_max_d, dp_max_b};
            return DB.at(ijk_idx).D;
        }
        if (ijk[2] == end) {
            DB[ijk_idx] = {inside({ijk[0], ijk[1]}) + outside(ijk) + inside({ijk[1], ijk[2]}), dp_max_b};
            return DB.at(ijk_idx).D;
        }
        std::array<int16_t, 3> cur_b;
        double cur_d;
        for (int16_t k = ijk[2] + 1; k < (int16_t)end+1; k++) {
            cur_b = {ijk[1], ijk[2], k};
            cur_d = inside({ijk[0], ijk[1]}) + outside(ijk) + dp(cur_b);
            if (cur_d > dp_max_d) {
                dp_max_d = cur_d;
                dp_max_b = cur_b;

            }
        }
        DB[ijk_idx] = {dp_max_d, dp_max_b};
        return DB.at(ijk_idx).D;
    };
    std::array<int16_t, 3> compare_lval_holder;
    double max_d = inside({(int16_t)start, (int16_t)end});
    std::array<int16_t, 3> max_b = {-1, -1, -1};
    for (int16_t j = start + 1; j < end; j++) {
        for (int16_t k = j + 1; k < end + 1; k++) {
            compare_lval_holder = {(int16_t)start, j, k};
            if (dp(compare_lval_holder) > max_d) {
                max_b = {(int16_t)start, j, k};
                max_d = dp(max_b);
            }
        }
    }

    return {DB, max_d, max_b, out_mem};
}



double get_high_threshold(uint32_t seg_len, const std::vector<double> & smoothed_threshold,
                          double threshold_rate)
{
    return (seg_len < smoothed_threshold.size()) ? smoothed_threshold[seg_len] : threshold_rate;
}

void refine_segmentation(const std::vector<double> &y_raw, std::vector<uint32_t> &y_idxs,
                                    double sigma, int skip, int min_interval_splice)
{
    std::vector<int> refine_y_idxs; refine_y_idxs.reserve(y_idxs.size());
    for (std::vector<uint32_t>::const_iterator it = y_idxs.begin(); it != --y_idxs.end();) {
        int s_yidx = *it;
        int e_yidx = *(++it);
        if (e_yidx - s_yidx <= 2*skip) continue;
        std::vector<double> i_vals; i_vals.reserve(e_yidx - s_yidx);
        for (int i = s_yidx; i < e_yidx; i++) i_vals.push_back(y_raw[i]);
        for (int i = 0; i < skip; i++) {
            i_vals[i] = 0.0f;
            i_vals[i_vals.size() - 1 - i] = 0.0f;
        }
        double sum = std::accumulate(i_vals.begin(), i_vals.end(), 0.0f);
        if ( sum < min_interval_splice) continue;
        std::vector<double> i_gauss = dirty_gaussian_filter1d(i_vals, sigma, 1.0f, true);
        std::deque<uint32_t> peaks = find_peaks(i_gauss.begin(), i_gauss.end(), i_gauss.begin(), i_gauss, 
                                                skip);
        for (auto & i : peaks) {
            int beg = (int)(i - sigma + 0.5f); if (beg < 0) beg = 0;
            int end = (int)(i + sigma + 1.5f); if (end > (int)i_gauss.size()) end = i_gauss.size();
            sum = std::accumulate(i_gauss.begin() + beg,
                            i_gauss.begin() + end, 0.0f);
            if (sum < min_interval_splice) continue;
            refine_y_idxs.push_back(i + s_yidx);
        }
    }
    y_idxs.insert(y_idxs.end(), refine_y_idxs.begin(), refine_y_idxs.end());
}

void get_unaligned_gaps_and_polyA(Read &read, const std::vector<std::array<int, 2>> &segs)
{

        
    read.gaps.reserve(200);
    bool one_in = false;
    for (auto & e : read.data) if (e == 1) {one_in = true; break;}
    if (!one_in) return;
    int d;
    std::vector<std::array<int, 2>> intervals; intervals.reserve(read.data.size());
    for (int i = 0; i < (int)read.data.size(); i++) {
        d = read.data[i];
        if (d != 1) continue;
        std::vector<std::array<int, 2>> group; group.reserve(read.data.size() - 1);
        while (d == 1 && i < (int)read.data.size()) {
            group.push_back({i, d});
            i++; d = read.data[i];
        }
        i--;
        int f_seg_idx = group[0][0];
        int l_seg_idx = group.back()[0];
        intervals.push_back({f_seg_idx, l_seg_idx});
    }

    int f_seg_idx = intervals[0][0];
    int start = segs[f_seg_idx][0];
    int l_seg_idx = intervals.back()[1];
    int end = segs[l_seg_idx][1];
    std::array<int, 2> q_ssc_pos = get_interval_start(start, read);
    std::array<int, 2> q_esc_pos = get_interval_end(end, read);


    int s, e, step;
    std::vector<Ilp> s_polys;
    for (auto & _char : char_holder) {
        s = 0;
        e = q_ssc_pos[0];
        step = 1;
        char sc_char = _char;
        std::vector<Ilp> ilp_queue;
        if (read.strand == '-') {
            s = read.seq.size() - s - 1;
            e = read.seq.size() - e - 1;
            step = -1;
            sc_char = getMap()[_char];
            ilp_queue = find_longest_poly_backwards(read.seq, s, e, step, sc_char);
        }
        else ilp_queue = find_longest_poly(read.seq, s, e, step, sc_char);
        s_polys.reserve(s_polys.size() + ilp_queue.size());
        for (auto & ilp : ilp_queue) {
            if (ilp.l < 20 || ilp.p < 0.85f) continue;
            s_polys.push_back(ilp);
            s_polys.back().c = _char;
        }

    }

    if (s_polys.size() > 0) {
        Ilp max_ilp; max_ilp.p = -std::numeric_limits<double>::max();
        for (auto & poly : s_polys) if (max_ilp.p < poly.p) max_ilp = poly;
        int poly_to_gene_gap_size = q_ssc_pos[0] - max_ilp.i - max_ilp.l;
        std::string temp_s = "S";
        temp_s += max_ilp.c;
        temp_s +=  "_" + std::to_string(max_ilp.l) + ":" + std::to_string(poly_to_gene_gap_size);
        read.gaps.emplace_back(temp_s);
        read.gaps.emplace_back("SSC:" + std::to_string(max_ilp.i));
    }
    else {
        read.gaps.emplace_back("SSC:" + std::to_string(q_ssc_pos[0]));
    }
    std::vector<Ilp> e_polys;
    for (auto & _char : char_holder) {
        s = q_esc_pos[0];
        e = read.length;
        step = 1;
        char sc_char = _char;
        std::vector<Ilp> ilp_queue;
        if (read.strand == '-') {
            s = read.seq.size() -s -1;
            e = read.seq.size() -e -1;
            step = -1;
            sc_char = getMap()[_char];
            ilp_queue = find_longest_poly_backwards(read.seq, s, e, step, sc_char);
        }
        else ilp_queue = find_longest_poly(read.seq, s, e, step, sc_char); 
        e_polys.reserve(e_polys.size() + ilp_queue.size());
        for (auto & ilp : ilp_queue) {
            if (ilp.l < 20 || ilp.p < 0.85f) continue;
            e_polys.push_back(ilp);
            if (read.strand == '-') {
                int afgagaobgag = 1;
            }
            e_polys.back().c = _char;
        }
        
    }
    if (e_polys.size() > 0) {

        Ilp max_ilp; max_ilp.p = -std::numeric_limits<double>::max();
        for (auto & poly : e_polys) if (max_ilp.p < poly.p) max_ilp = poly;
        int poly_to_gene_gap_size = max_ilp.i;
        std::string temp_s = "E";
        temp_s += max_ilp.c;
        temp_s +=  "_" + std::to_string(max_ilp.l) + ":" + std::to_string(poly_to_gene_gap_size);
        read.gaps.emplace_back(std::move(temp_s));
        read.gaps.emplace_back("ESC:" + std::to_string(read.length - q_esc_pos[0] - poly_to_gene_gap_size));
    }
    else {
        read.gaps.emplace_back("ESC:" + std::to_string(read.length - q_esc_pos[0]));
    }
    for (auto it = intervals.begin(); it != intervals.end() -  1;) {
        int il_l_seg_idx = (*it)[1];
        int il_end = segs[il_l_seg_idx][1];
        std::array<int, 2> q_gap_s_and_s_slack = get_interval_end(il_end, read);
        int i2_f_seg_idx = (*(++it))[0];
        int i2_start = segs[i2_f_seg_idx][0];
        std::array<int, 2> q_gap_e_and_e_slack = get_interval_start(i2_start, read);
        int q_gap_size = q_gap_e_and_e_slack[0] - q_gap_s_and_s_slack[0];
        q_gap_size = std::max(0, q_gap_size + q_gap_e_and_e_slack[1] + q_gap_s_and_s_slack[1]);
        read.gaps.emplace_back(std::to_string(il_l_seg_idx) + "-" + std::to_string(i2_f_seg_idx) + ":" +
                         std::to_string(q_gap_size));
    }
    fast_set(read.gaps);
}

std::array<int, 2> get_interval_start(int start, const Read &read)
{
    for (const auto & e : read.intervals) {
        int t_start = e.first[0];
        int t_end = e.first[1];
        int q_start = e.first[2];
        int q_end = e.first[3];
        const auto & cigar = e.second;
        if (t_end < start) continue;
        int q_pos, slack;
        if (start < t_start) {
            q_pos = q_start;
            slack = start - t_start;
        }
        else {
            q_pos = forward_thread_cigar(cigar, start, t_start, q_start);
            slack = 0;
        }
        return {q_pos, slack};
    }
    assert(false);
    return {};
}

int forward_thread_cigar(const std::vector<cig_element> &cigar, int t_goal, int t_pos, int q_pos)
{
    int cig_idx = 0; int c; char t;
    while (t_pos < t_goal) {
        t = cigar[cig_idx].c; c = cigar[cig_idx].count;
        c = std::min(c, t_goal - t_pos);
        switch (t) {
        case 'M':
        case 'X':
        case '=':
            t_pos += c;
            q_pos += c;
            break;
        case 'D':
            t_pos += c;
            break;
        case 'I':
            q_pos += c;
            break;
        }
        cig_idx += 1;
    }
    assert(t_pos == t_goal);
    return q_pos;
}

std::vector<Ilp> find_longest_poly_backwards(const std::string &seq, int s, int e, int step,
                                                 char _char, int match_score, int mismatch_score)
{
    if (e - s == 0) return {};
    std::vector<int> scores(1); scores.reserve((s - e) / (-1*step) + 1);
    if (seq[s] == _char) scores[0] = match_score;
    else scores[0] = 0;
    for (int i = s + step; i > e; i+=step) {
        scores.push_back(std::max(0, scores.back() + ((seq[i] == _char) ? match_score : mismatch_score)));
    }
    bool k;
    std::vector<Ilp> stack; stack.reserve(scores.size());
    for (int i = 0; i < (int)scores.size(); i++) {
        k = scores[i] > 0; 
        if (!k) continue;
        std::deque<std::array<int, 2>> group;
        while (k == 1) {
            group.push_back({i, scores[i]});
            i++;  
            if (i < (int)scores.size()) k = scores[i] > 0;
            else k = 0;
        }
        int max_i, max_s; max_s = 0; max_i = 0;
        for (auto & g : group) {
            if (g[1] > max_s) {
                max_s = g[1];
                max_i = g[0];
            }
        }
       int l = max_i + 1 - group[0][0];
       std::string temp_s; temp_s.reserve((s - e) / (-1*step) + 1);
       for (int j = s; j > e; j += step) temp_s.push_back(seq[j]);
       float count = 0;
       
       for (int j = group[0][0]; j < group[0][0] + l; j++) if (temp_s[j] == _char) count += 1.0f;
       stack.push_back({group[0][0], l, count/l, 0});
    }
    stack.shrink_to_fit();
    return stack;
}


std::vector<Ilp> find_longest_poly(const std::string &seq, int s, int e, int step,
                                                 char _char, int match_score, int mismatch_score)
{
    if (e - s == 0) return {};
    std::vector<int> scores(1); scores.reserve((e - s) / step + 1);
    if (seq[s] == _char) scores[0] = match_score;
    else scores[0] = 0;
    for (int i = s + step; i < e; i+=step) {
        scores.push_back(std::max(0, scores.back() + ((seq[i] == _char) ? match_score : mismatch_score)));
    }
    bool k;
    std::vector<Ilp> stack; stack.reserve(scores.size());
    for (int i = 0; i < (int)scores.size(); i++) {
        k = scores[i] > 0; 
        if (!k) continue;
        std::deque<std::array<int, 2>> group;
        while (k == 1) {
            group.push_back({i, scores[i]});
            i++;  
            if (i < (int)scores.size()) k = scores[i] > 0;
            else k = 0;
        }
        int max_i, max_s; max_s = 0; max_i = 0;
        for (auto & g : group) {
            if (g[1] > max_s) {
                max_s = g[1];
                max_i = g[0];
            }
        }
       int l = max_i + 1 - group[0][0];
       std::string temp_s; temp_s.reserve((e - s) / step + 1);
       for (int j = s; j < e; j +=step) temp_s.push_back(seq[j]);
       float count = 0;
       
       for (int j = group[0][0]; j < group[0][0] + l; j++) if (temp_s[j] == _char) count += 1.0f;
       stack.push_back({group[0][0], l, count/l, 0});
    }
    stack.shrink_to_fit();
    return stack;
}

std::array<int, 2> get_interval_end(int end, const Read &read)
{
    for (auto e = read.intervals.rbegin(); e != read.intervals.rend(); ++e) {
        int t_start = e->first[0];
        int t_end = e->first[1];
        int q_start = e->first[2];
        int q_end = e->first[3];
        const auto & cigar = e->second;
        if (t_start > end) continue;
        int q_pos, slack;
        if (t_end < end) {
            q_pos = q_end;
            slack = t_end - end;
        }
        else {
            q_pos = forward_thread_cigar(cigar, end, t_start, q_start);
            slack = 0;
        }
        return {q_pos, slack};
    }
    assert(false);
    return {};
}
