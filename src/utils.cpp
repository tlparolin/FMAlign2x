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

// FMAlign2x - An extended version of FMAlign2 for aligning multiple ultra-long sequences
// Author: Thiago Luiz Parolin
// Contact: thiago.parolin@unesp.br
// July 2025

#include "utils.h"

KSEQ_INIT(int, read)

/**
 * @brief A timer class that measures elapsed time.
 * This class uses C++11 chrono library to measure elapsed time in seconds with double precision.
 * The timer starts at construction and can be reset to zero by calling reset().
 * The elapsed time can be obtained by calling elapsed_time() method.
 * The timer is based on std::chrono::steady_clock, which is a monotonic clock that is not subject to system clock adjustments.
 */
Timer::Timer() noexcept : start_time_{std::chrono::steady_clock::now()} {}

void Timer::reset() noexcept { start_time_ = std::chrono::steady_clock::now(); }

[[nodiscard]] double Timer::elapsed_time() const noexcept {
    return std::chrono::duration<double>(std::chrono::steady_clock::now() - start_time_).count();
}

/**
 * @brief: read fasta and fastq format data
 * @param data_path   the path to the target data
 * @param data store sequence content
 * @param name store sequence name
 * @return multiple sequence stored in vector
 */
void read_data(const std::filesystem::path &data_path, std::vector<std::string> &data, std::vector<std::string> &name, bool verbose) {
    if (verbose && global_args.verbose) {
        std::cout << "#                   Reading Data...                         #" << std::endl;
        print_table_divider();
    }

    if (!access_file(data_path)) {
        print_table_bound();
        std::cerr << "Error: " << data_path.string() << " could not be accessed. Please check the path or if the data exists!" << std::endl
                  << "Program Exit!" << std::endl;
        std::exit(1);
    }

    if (verbose && global_args.verbose) {
        print_table_line(std::format("{} could be accessed", data_path.string()));
    }

    FILE *f_pointer = std::fopen(data_path.string().c_str(), "r");
    if (!f_pointer) {
        throw std::runtime_error(std::format("Unable to open file: {}", data_path.string()));
    }
    kseq_t *file_t = kseq_init(fileno(f_pointer));

    uint64_t merged_length = 0;
    int64_t tmp_length = 0;
    // stop loop when tmp_length equals -1
    while ((tmp_length = kseq_read(file_t)) >= 0) {
        std::string tmp_data = clean_sequence(file_t->seq.s);
        std::string tmp_name = file_t->name.s;
        if (file_t->comment.s)
            tmp_name += file_t->comment.s;

        data.emplace_back(std::move(tmp_data));
        name.emplace_back(std::move(tmp_name));
        merged_length += tmp_length;
    }

    kseq_destroy(file_t);
    std::fclose(f_pointer);

    if (verbose && global_args.verbose && merged_length + data.size() > UINT32_MAX && M64 == 0) {
        print_table_bound();
        std::cerr << "Error: The input data is too large and the 32-bit program may not produce correct results." << std::endl
                  << "Please compile a 64-bit program using the M64 parameter." << std::endl
                  << "Program Exit!" << std::endl;
        std::exit(1);
    }

#if M64
    if (verbose && global_args.verbose) {
        print_table_line(std::format("Data Memory Usage: {:.2f} GB", merged_length / double(1ULL << 30)));
    }
#else
    if (verbose && global_args.verbose) {
        print_table_line(std::format("Data Memory Usage: {:.2f} MB", merged_length / double(1ULL << 20)));
    }
#endif
    if (verbose && global_args.verbose) {
        print_table_line(std::format("Sequence Number: {}", data.size()));
        print_table_divider();
    }
}

/**
 * @brief: Check whether the file exists in the specified path.
 * @param data_path   The file path to check.
 * @return Returns true if the file exists, otherwise false.
 */
bool access_file(const std::filesystem::path &file_path) {
    std::ifstream file(file_path);
    return file.good();
}

/**
 * @brief Cleans the input sequence by removing any non-ATCG characters.
 * This function removes any characters from the input sequence that are not A, T, C, or G (case-insensitive).
 * The cleaned sequence is then returned as a new string.
 * @param sequence The input DNA sequence to be cleaned.
 * @return The cleaned DNA sequence as a new string.
 */
std::string clean_sequence(std::string sequence) {
    std::string result(sequence.size(), '-');
    std::transform(std::execution::par, sequence.begin(), sequence.end(), result.begin(), [](char c) -> char {
        switch (c) {
        case 'a':
        case 'A':
            return 'A';
        case 'c':
        case 'C':
            return 'C';
        case 'g':
        case 'G':
            return 'G';
        case 't':
        case 'T':
            return 'T';
        case 'u':
        case 'U':
            return 'U';
        default:
            return '-';
        }
    });
    return result;
}

void ArgParser::add_argument(const std::string &name, bool required, const std::string &default_value) {
    if (args_.contains(name)) {
        throw std::invalid_argument("Duplicate argument name: " + name);
    }
    args_[name] = {required, default_value};
}

void ArgParser::add_argument_help(const std::string &name, const std::string &help_text) {
    if (!args_.contains(name)) {
        throw std::invalid_argument("Invalid argument name: " + name);
    }
    args_[name].help_text = help_text;
}

void ArgParser::print_help() const {
    std::cout << "Usage: FMAlign2 [OPTIONS]" << std::endl << std::endl << "Options:" << std::endl;
    for (const auto &[key, arg] : args_) {
        std::cout << "  -" << key;
        if (!arg.required) {
            std::cout << " (optional)";
        }
        if (!arg.default_value.empty()) {
            std::cout << " [default: " << arg.default_value << "]";
        }
        std::cout << "\n    " << arg.help_text << "\n\n";
    }
}

void ArgParser::parse_args(int argc, char **argv) {
    for (int i = 1; i < argc; i++) {
        std::string arg = argv[i];
        if (arg.starts_with('-')) {
            std::string name = arg.substr(1);
            if (name == "help" || name == "h") {
                print_help();
                std::exit(0);
            }
            if (!args_.contains(name)) {
                throw std::invalid_argument("Invalid argument: " + arg);
            }
            Arg &a = args_[name];
            if (!a.value.empty()) {
                throw std::invalid_argument("Duplicate argument: " + arg);
            }
            if (i == argc - 1 || argv[i + 1][0] == '-') {
                if (a.required) {
                    throw std::invalid_argument("Missing value for argument: " + arg);
                }
                a.value = a.default_value;
            } else {
                a.value = argv[++i];
            }
        }
    }
    for (auto &[key, arg] : args_) {
        if (arg.required && arg.value.empty()) {
            throw std::invalid_argument("Missing required argument: -" + key);
        }
        if (!arg.required && arg.value.empty()) {
            arg.value = arg.default_value;
        }
    }
}

std::string ArgParser::get(const std::string &name) const {
    if (!args_.contains(name)) {
        throw std::invalid_argument("Invalid argument name: " + name);
    }
    const Arg &a = args_.at(name);
    if (a.value.empty()) {
        throw std::invalid_argument("Missing value for argument: --" + name);
    }
    return a.value;
}

bool ArgParser::has(const std::string &name) const { return args_.contains(name) && !args_.at(name).value.empty(); }

/**
 * @brief Print information about the FMAlign2 algorithm
 * This function prints various information about the FMAlign2 algorithm,
 * including the mode (32-bit or 64-bit), number of threads, minimum MEM length,
 * sequence coverage, and parallel align method.
 * @return void
 */
void print_algorithm_info() {
    print_table_bound();
    std::cout << "#               FMAlign2x algorithm info                    #\n";
    print_table_divider();

#if M64
    print_table_line("Mode: 64 bit");
#else
    print_table_line("Mode: 32 bit");
#endif
    print_table_line(std::format("Thread: {}", global_args.thread));

    if (global_args.min_mem_length < 0) {
        print_table_line("Minimum MEM length: square root of mean length");
    } else {
        print_table_line(std::format("Minimum MEM length: {}", global_args.min_mem_length));
    }

    if (global_args.min_seq_coverage < 0) {
        print_table_line("Sequence coverage: default");
    } else {
        print_table_line(std::format("Sequence coverage: {}", global_args.min_mem_length));
    }

    print_table_line(std::format("Parallel align method: {}", global_args.package));

    print_table_bound();
}

void print_table_line(const std::string &output) {
    std::cout << "# " << std::left << std::setw(TABLE_LEN - 2) << std::setfill(' ') << output << "#" << std::endl;
}

void print_table_divider() { std::cout << "#" << std::right << std::setw(TABLE_LEN) << std::setfill('-') << "#" << std::endl; }

void print_table_bound() { std::cout << "#" << std::right << std::setw(TABLE_LEN) << std::setfill('#') << "#" << std::endl; }
