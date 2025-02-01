/*
 * Copyright [2023] [MALABZ_UESTC Pinglu Zhang]
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

// Author: Pinglu Zhang
// Contact: zpl010720@gmail.com
// Created: 2023-02-25

// The main function is the entry point of the program. It is where the program starts executing. 
// the program starts executing. 
#include "include/common.h"
#include "include/utils.h"
#include "include/mem_finder.h"
#include "include/sequence_split_align.h"
#if defined(__linux__)
#include "include/thread_pool.h"
#endif
#include <sharg/all.hpp>
#include <thread>
#include <stdexcept>

GlobalArgs global_args;

void initialise_parser(sharg::parser & parser)
{
    parser.info.author = "FMAlign2 SeqAn version - Thiago Luiz Parolin";
    parser.info.short_description = "FMAlign2 using SeqAn library is based on original work from Pinglu Zhang";
    parser.info.version = "0.0.1";
    parser.info.date = "31-01-2025";
}

int main(int argc, char** argv) {
    Timer timer;
    sharg::parser parser{"FMAlign2", argc, argv};
    initialise_parser(parser);

    global_args.data_path = "data/mt1x.fasta";
    global_args.thread = std::thread::hardware_concurrency();
    global_args.min_mem_length = -1;
    global_args.min_seq_coverage = -1;
    sharg::arithmetic_range_validator seq_coverage_validator{0.1, 1.0};
    global_args.package = "mafft";
    sharg::value_list_validator package_validator{"mafft", "halign2", "halign3"};
    global_args.output_path = "output.fmaligned2.fasta";
    global_args.degree = 0;
    global_args.filter_mode = "default";
    sharg::value_list_validator filter_validator{"global", "local", "default"};
    global_args.verbose = 1;
    sharg::value_list_validator verbose_validator{0, 1};
    
    parser.add_option(global_args.data_path, sharg::config{.short_id = 'i', .long_id = "input", .description = "The path to the input file.", .required = true});
    parser.add_option(global_args.thread, sharg::config{.short_id = 't', .long_id = "threads", .description = "Max number of threads (default: number of CPUs)."});
    parser.add_option(global_args.min_mem_length, sharg::config{.short_id = 'l', .long_id = "min-mem-length", .description = "Minimum length of MEM (default: sqrt of mean length)."});
    parser.add_option(global_args.min_seq_coverage, sharg::config{.short_id = 'c', .long_id = "coverage", .description = "A floating-point parameter that specifies the minimum coverage across all sequences, with values ranging from 0 to 1. The default \
setting is that if sequence number less 100, parameter is set to 1 otherwise 0.7.", .validator = seq_coverage_validator});
    parser.add_option(global_args.package, sharg::config{.short_id = 'p', .long_id = "package", .description = "MSA method: halign2, halign3, or mafft.", .validator = package_validator});
    parser.add_option(global_args.output_path, sharg::config{.short_id = 'o', .long_id = "output", .description = "Path to the output file."});
    parser.add_option(global_args.degree, sharg::config{.short_id = 'd', .long_id = "depth", .description = "Depth of recursion (ignored)."});
    parser.add_option(global_args.filter_mode, sharg::config{.short_id = 'f', .long_id = "filter-mode", .description = "Filter MEMs mode: global or local (default: auto).", .validator = filter_validator});
    parser.add_option(global_args.verbose, sharg::config{.short_id = 'v', .long_id = "verbose", .description = "Verbose option (0 or 1).", .validator = verbose_validator});
    
    try {
        parser.parse();

    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\nProgram Exit!" << std::endl;
        return 1;
    }

    if (global_args.verbose) {
        print_algorithm_info();
    }

    std::vector<std::vector<seqan3::dna4>> data;

    try {
        read_data(global_args.data_path.c_str(), data, true);
        // std::vector<std::vector<std::pair<int_t, int_t>>> split_points_on_sequence = find_mem(data);
        // split_and_parallel_align(data, name, split_points_on_sequence);
    } catch (const std::bad_alloc& e) {
        print_table_bound();
        std::cerr << "Error: " << e.what() << "\nProgram Exit!" << std::endl;
        std::exit(1);
    }

    double total_time = timer.elapsed_time();
    if (global_args.verbose) {
        std::cout << "FMAlign2 total time: " << std::fixed << std::setprecision(2) << total_time << " seconds." << std::endl;
    }

    return 0;
}
