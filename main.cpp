

#include "argparser.h"
#include "readSplitcsv.h"
#include <chrono>
#include <atomic>
#include <iostream>
#include <deque>
#include <filesystem>
#include "freddie_segment.h"
#include <future>
#include <deque>
namespace fs = std::filesystem;
int main(int argc, char** argv)
{
    ArgParser arg = arg_parse(argc, argv);
    fs::create_directory(arg.out_dir);
    std::deque<segArgs> segment_args;

    for (const auto & contig : fs::directory_iterator(arg.split_dir)) {
        if (!contig.is_directory()) continue;
        fs::create_directory(arg.out_dir + "/" + contig.path().filename().string());
        for (const auto & file : fs::directory_iterator(contig)) {
            const std::string & FilenameAlias = file.path().filename().string();
            if (!(FilenameAlias.substr(0,6) == "split_" &&
                    FilenameAlias.substr(FilenameAlias.size() - 4, 4) == ".tsv")) continue;

            int idx = FilenameAlias.find_last_of("_") + 1;
            int o_idx = FilenameAlias.find_last_of(".") - idx;
            uint32_t tint_id = std::stoi(FilenameAlias.substr(idx, o_idx));

            segment_args.emplace_back(segArgs(arg, smooth_threshold(arg.threshold_rate),
                                              contig.path().filename().string(), tint_id));

        }
    }
    uint32_t idx = 0;

    if (arg.threads == 1) {

        for (const auto& seg : segment_args) {
            run_segment(0, seg);
            if ((idx % ((int)(segment_args.size() + 100 - 1) / 100 ))) {
                idx++;
                continue;
            }
        
            printf("[freddie_segment] Done with %u/%u tints (%.1f)%% \n", idx, (uint32_t)segment_args.size(),
                (float)idx / segment_args.size() * 100);
            idx++;
        }
    }
    else {
        ctpl::thread_pool p(arg.threads);
        int thread_counter = 0;
        std::deque<std::future<void>> futures_queue;

        for (const auto& seg : segment_args) {
            futures_queue.push_back(p.push(run_segment, seg)); 
        }
        for (auto & future : futures_queue) {
            future.get();
            if ((idx % ((int)(segment_args.size() + 100 - 1) / 100 ))) {
                idx++;
                continue;
            }
        
            printf("[freddie_segment] Done with %u/%u tints (%.1f)%% \n", idx, (uint32_t)segment_args.size(),
                (float)idx / segment_args.size() * 100);
            idx++;
        }
    }


    return 0;
}
