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
#include <thread>
#include "mpi.h"

GlobalArgs global_args;
int main(int argc, char** argv) {

   // MPI initialization
    MPI_Init(&argc, &argv);

    int world_rank, world_size;

    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    // Create a Timer object to record the execution time.
    Timer timer;

    // Create an ArgParser object to parse command line arguments.
    ArgParser parser;
    std::string output = "";

    // Add command line arguments to the ArgParser object.
    parser.add_argument("i", true, "data/mt1x.fasta");
    parser.add_argument_help("i", "The path to the input file.");
    parser.add_argument("t", false, "cpu_num");
    parser.add_argument_help("t", "The maximum number of threads that the program runs, the recommended setting is the number of CPUs.");
    parser.add_argument("l", false, "default");
    parser.add_argument_help("l", "The minimum length of MEM, the default value is square root of mean length.");
    parser.add_argument("c", false, "1");
    parser.add_argument_help("c", "A floating-point parameter that specifies the minimum coverage across all sequences, with values ranging from 0 to 1. The default \
setting is that if sequence number less 100, parameter is set to 1 otherwise 0.7.");
    parser.add_argument("p", false, "mafft");
    parser.add_argument_help("p", "The MSA method used in parallel align. for example, halign3, halign2 and mafft.");
    parser.add_argument("o", false, "output.fmaligned2.fasta");
    parser.add_argument_help("o", "The path to the output file.");
    parser.add_argument("d", false, "0");
    parser.add_argument_help("d", "Depth of recursion, you could ignore it.");
    parser.add_argument("f", false, "default");
    parser.add_argument_help("f", "The filter MEMs mode. The default setting is that if sequence number less 100, local mode otherwise global mode.");
    parser.add_argument("v", false, "1");
    parser.add_argument_help("v", "Verbose option, 0 or 1. You could ignore it.");
    parser.add_argument("h", false, "help");
    parser.add_argument_help("h", "print help information");

    // Add command line arguments to the ArgParser object.
    try {
        parser.parse_args(argc, argv);
        global_args.data_path = parser.get("i");
        std::string tmp_thread = parser.get("t");
        if (tmp_thread == "cpu_num") {
            global_args.thread = world_size;
        }
        else {
            global_args.thread = std::stoi(tmp_thread);
        }
        std::string tmp_len = parser.get("l");
        if (tmp_len != "default") {
            global_args.min_mem_length = std::stoi(parser.get("l"));
        }
        else {
            global_args.min_mem_length = -1;
        }

        std::string tmp_filter_mode = parser.get("f");
        if (tmp_filter_mode == "default") {
            global_args.filter_mode = tmp_filter_mode;
        }
        else if (tmp_filter_mode == "global" || tmp_filter_mode == "local") {
            global_args.filter_mode = tmp_filter_mode;
        }
        else {
            throw "filer mode --f parameter should be global or local!";
        }

        global_args.verbose = std::stoi(parser.get("v"));
        if (global_args.verbose != 0 && global_args.verbose != 1) {
            throw "verbose should be 1 or 0";
        }

        global_args.degree = std::stoi(parser.get("d"));
        if (global_args.degree > 2) {
            exit(1);
        }

        std::string tmp_c = parser.get("c");
        if (tmp_c == "default") {
            global_args.min_seq_coverage = -1;
        }
        else {
            global_args.min_seq_coverage = std::stof(parser.get("c"));
            if (global_args.min_seq_coverage < 0 || global_args.min_seq_coverage > 1) {
                throw "Error: min_seq_coverage should be ranged from 0 to 1!";
            }
        }

        global_args.package = parser.get("p");
        if (global_args.package != "halign2" && global_args.package != "halign3" && global_args.package != "mafft") {
            throw ("Error: " + global_args.package + " is a invalid method!");
        }

        global_args.output_path = parser.get("o");
    } // Catch any invalid arguments and print the help message.
    catch (const std::invalid_argument& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        std::cerr << "Program Exit!" << std::endl;
        parser.print_help();
        return 1;
    }

    if (world_rank == 0) {
        if (global_args.verbose) {
            print_algorithm_info(world_size);
        }
    }

    std::vector<std::string> data;
    std::vector<std::string> name;
    std::vector<std::vector<std::pair<int_t, int_t>>> split_points_on_sequence;

    try {
        if (world_rank == 0){
            // Read data from the input file and store in data and name vectors
            read_data(global_args.data_path.c_str(), data, name, true);

            // Find MEMs in the sequences and split the sequences into fragments for parallel alignment.
            split_points_on_sequence = find_mem(data);
        }

        // Broadcast data size
        size_t data_size = data.size();
        MPI_Bcast(&data_size, 1, MPI_UNSIGNED_LONG, 0, MPI_COMM_WORLD);

        // Resize data vector for other ranks
        if (world_rank != 0) {
            data.resize(data_size);
        }

        // Broadcast each string in the data vector
        for (size_t i = 0; i < data_size; ++i) {
            size_t str_size = data[i].size();
            MPI_Bcast(&str_size, 1, MPI_UNSIGNED_LONG, 0, MPI_COMM_WORLD);

            if (world_rank != 0) {
                data[i].resize(str_size);
            }

            MPI_Bcast(data[i].data(), str_size, MPI_CHAR, 0, MPI_COMM_WORLD);
        }

        // Broadcast name size
        size_t name_size = name.size();
        MPI_Bcast(&name_size, 1, MPI_UNSIGNED_LONG, 0, MPI_COMM_WORLD);

        // Resize name vector for other ranks
        if (world_rank != 0) {
            name.resize(name_size);
        }

        // Broadcast name strings
        for (size_t i = 0; i < name_size; ++i) {
            size_t str_size = name[i].size();
            MPI_Bcast(&str_size, 1, MPI_UNSIGNED_LONG, 0, MPI_COMM_WORLD);

            if (world_rank != 0) {
                name[i].resize(str_size);
            }

            MPI_Bcast(name[i].data(), str_size, MPI_CHAR, 0, MPI_COMM_WORLD);
        }

        // Broadcast split_points_on_sequence size
        size_t split_size = split_points_on_sequence.size();
        MPI_Bcast(&split_size, 1, MPI_UNSIGNED_LONG, 0, MPI_COMM_WORLD);

        // Resize split_points_on_sequence for other ranks
        if (world_rank != 0) {
            split_points_on_sequence.resize(split_size);
        }

        for (size_t i = 0; i < split_size; ++i) {
            size_t pair_count = split_points_on_sequence[i].size();
            MPI_Bcast(&pair_count, 1, MPI_UNSIGNED_LONG, 0, MPI_COMM_WORLD);

            if (world_rank != 0) {
                split_points_on_sequence[i].resize(pair_count);
            }

            MPI_Bcast(split_points_on_sequence[i].data(), pair_count * sizeof(std::pair<int_t, int_t>), MPI_BYTE, 0, MPI_COMM_WORLD);
        }

        // All nodes execute split_and_parallel_align
        split_and_parallel_align(data, name, split_points_on_sequence, world_rank, world_size); 
    }
    catch (const std::bad_alloc& e) { // Catch any bad allocations and print an error message.
        print_table_bound();
        std::cerr << "Error: " << e.what() << std::endl;
        std::cout << "Program Exit!" << std::endl;
        exit(1);
    }

        double total_time = timer.elapsed_time();
        std::stringstream s;
        s << std::fixed << std::setprecision(2) << total_time;
        if (global_args.verbose) {
            output = "FMAlign2 total time: " + s.str() + " seconds.";
            print_table_line(output);
            print_table_bound();
        }

    // Finalizar MPI
    MPI_Finalize();

    return 0;
}
