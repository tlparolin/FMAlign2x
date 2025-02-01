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
// Created: 2023-02-24


#include "../include/utils.h"

/**
 * @brief A timer class that measures elapsed time. 
 * This class uses C++11 chrono library to measure elapsed time in seconds with double precision. 
 * The timer starts at construction and can be reset to zero by calling reset(). 
 * The elapsed time can be obtained by calling elapsed_time() method. 
 * The timer is based on std::chrono::steady_clock, which is a monotonic clock that is not subject to system clock adjustments.
*/
Timer::Timer() {
    start_time_ = std::chrono::steady_clock::now();
}
void Timer::reset() {
    start_time_ = std::chrono::steady_clock::now();
}

double Timer::elapsed_time() const {
    std::chrono::duration<double> elapsed = std::chrono::steady_clock::now() - start_time_;
    return elapsed.count();
}

/**
 * @brief: read fasta and fastq format data
 * @param data_path   the path to the target data
 * @param data store sequence content
 * @param name store sequence name
 * @return multiple sequence stored in vector 
*/
void read_data(const char* data_path, std::vector<std::vector<seqan3::dna4>>& data, bool verbose){
    if (verbose && global_args.verbose) {
        std::cout << "#                   Reading Data...                         #" << std::endl;
        print_table_divider();
    }
    std::string output = "";
    std::string str_data_path = data_path;

    // Load file using _dna4 alphabet
    seqan3::sequence_file_input file_in{str_data_path};
    uint64_t merged_length = 0;
    size_t sequence_count = 0;
    for (auto & record : file_in) {
        sequence_count++;
        merged_length += record.sequence().size();
    }

    if(verbose&& global_args.verbose && merged_length + sequence_count > UINT32_MAX && M64 == 0){
        print_table_bound();
        std::cerr << "Error: The input data is too large and the 32-bit program may not produce correct results. Please compile a 64-bit program using the M64 parameter." << std::endl;
        std::cerr << "Program Exit!" << std::endl;
        exit(1);
    }
#if M64
    if (verbose && global_args.verbose) {
        std::stringstream s;
        s << std::fixed << std::setprecision(2) << merged_length / pow(2, 30);
        output = "Data Memory Usage: " + s.str() + " GB";
        print_table_line(output);
    }
    
#else
    if (verbose && global_args.verbose) {
        std::stringstream s;
        s << std::fixed << std::setprecision(2) << merged_length / pow(2, 20);
        output = "Data Memory Usage: " + s.str() + " MB";
        print_table_line(output);
    }
#endif
    if (verbose && global_args.verbose) {
        std::string output = "Sequence Number: " + std::to_string(sequence_count);
        print_table_line(output);
        print_table_divider();
    }
    return;
}

/**
* @brief Print information about the FMAlign2 algorithm
* This function prints various information about the FMAlign2 algorithm,
* including the mode (32-bit or 64-bit), number of threads, minimum MEM length,
* sequence coverage, and parallel align method.
* @return void
*/
void print_algorithm_info() {
    print_table_bound();
    std::cout << "#               FMAlign2 algorithm info                     #" << std::endl;
    print_table_divider();
#if M64
    std::string output = "Mode: 64 bit";
    print_table_line(output);
#else
    std::string output = "Mode: 32 bit";
    print_table_line(output);
#endif
    std::string thread_output = "Thread: " + std::to_string(global_args.thread);
    print_table_line(thread_output);

    std::string l_output;
    if (global_args.min_mem_length < 0) {
        l_output = "Minimum MEM length: square root of mean length";
    }
    else {
        l_output = "Minimum MEM length: " + std::to_string(global_args.min_mem_length);
    }
    
    print_table_line(l_output);

    std::stringstream s;
    s << std::fixed << std::setprecision(2) << global_args.min_seq_coverage;

    std::string c_output;
    if (global_args.min_seq_coverage < 0) {
        c_output = "Sequence coverage: default";
    }
    else {
        c_output = "Sequence coverage: " + std::to_string(global_args.min_mem_length);
    }

    print_table_line(c_output);

    std::string p_output = "Parallel align method: " + global_args.package;
    print_table_line(p_output);

    print_table_bound();
}

void print_table_line(const std::string &output) {
    std::cout << "# " << std::left << std::setw(TABLE_LEN-2) << std::setfill(' ') << output << "#" << std::endl;
}

void print_table_divider() {
    std::cout << "#" << std::right << std::setw(TABLE_LEN) << std::setfill('-') << "#" << std::endl;
}

void print_table_bound() {
    std::cout << "#" << std::right << std::setw(TABLE_LEN) << std::setfill('#') << "#" << std::endl;
}


