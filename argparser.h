#pragma once
#include <string>
#include <vector>
#include <deque>
#include <algorithm>
#include <unordered_map>
#include <sstream>
#include <filesystem>

static const std::unordered_map<std::string, int> FLAGS = {{"-s", 0}, {"--split-dir", 0},
                                                    {"-o", 1}, {"--outdir", 1},
                                                    {"-t", 2}, {"--threads", 2},
                                                    {"-sd", 3}, {"--sigma", 3},
                                                    {"-tp", 4}, {"--threshold-rate", 4},
                                                    {"-vf", 5}, {"--variance-factor", 5},
                                                    {"-mps", 6}, {"--max-problem-size", 6},
                                                    {"-lo", 7}, {"--min-read-support-outside", 7}
                                              };


struct ArgParser {
    std::string split_dir;
    std::string out_dir;
    int threads;
    double sigma;
    double threshold_rate;
    double variance_factor;
    uint32_t max_problem_size;
    int min_read_support_outside;
    ArgParser() : out_dir("../freddie_segment"), threads(1), sigma(5.0f), threshold_rate(0.9f), variance_factor(3.0f),
        max_problem_size(50), min_read_support_outside(3) {}
};



struct segArgs : ArgParser {
    std::vector<double> smooth_threshold;
    std::string contig;
    uint32_t tint_id;
    segArgs(const ArgParser & baseArgs, const std::vector<double> && s_t,
            const std::string& con, const uint32_t tint)
        : ArgParser(baseArgs), smooth_threshold(s_t), contig(con), tint_id(tint) {}
};

inline ArgParser arg_parse(int argc, char ** argv_p) {
    ArgParser arg;
//    argc = 7;
   std::vector<std::string> argv(argv_p, argv_p + argc); 
//    argv[1] = std::string("-s");
//    argv[2] = "test/split/"; argv[3] =  "--outdir"; argv[4] = "test/freddie_segment/"; argv[5] =  "-t"; argv[6] =  "1";
    for (int i = 1; i < argc; i+=2) {
        int idx = FLAGS.at(argv[i]);
        switch(idx) {
        case 0:
            arg.split_dir = std::filesystem::current_path().string() + "/" + argv[i+1];
            break;

        case 1:
            arg.out_dir = std::filesystem::current_path().string() + "/" + argv[i+1];
            break;

        case 2:
            arg.threads = std::stoi(argv[i+1]);
            break;

        case 3:
            arg.sigma = std::stod(argv[i+1]);
            break;

        case 4:
            arg.threshold_rate= std::stod(argv[i+1]);
            break;

        case 5:
            arg.variance_factor = std::stod(argv[i+1]);
            break;

        case 6:
            arg.max_problem_size = std::stoi(argv[i+1]);
            break;

        case 7:
            arg.min_read_support_outside = std::stoi(argv[i+1]);
            break;
        }
    }
    return arg;
}


