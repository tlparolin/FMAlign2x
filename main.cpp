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

    if (world_rank == 0 && global_args.verbose) {
        print_algorithm_info(world_size, world_rank);
    }

    std::vector<std::string> data;
    std::vector<std::string> name;

    try {
        // Read data from the input file and store in data and name vectors
        read_data_mpi(global_args.data_path.c_str(), data, name, world_rank, world_size, global_args.verbose);

        // Find MEMs in the sequences and split the sequences into fragments for parallel alignment.
        std::vector<std::vector<std::pair<int_t, int_t>>> split_points_on_sequence;
        if (world_rank == 0) {
            split_points_on_sequence = find_mem(data, world_rank, world_size);
        }

        // Flatten the split_points_on_sequence into a continuous array for broadcasting
        std::vector<std::pair<int_t, int_t>> flattened_data;
        std::vector<size_t> sizes;
        if (world_rank == 0) {
            // Flatten the data into a single vector
            size_t total_size = 0;
            for (const auto& vec : split_points_on_sequence) {
                total_size += vec.size();
            }
            flattened_data.resize(total_size);

            size_t idx = 0;
            for (const auto& vec : split_points_on_sequence) {
                for (const auto& elem : vec) {
                    flattened_data[idx++] = elem;
                }
            }

            // Store the sizes of each subvector for later unflattening
            sizes.resize(split_points_on_sequence.size());
            for (size_t i = 0; i < split_points_on_sequence.size(); i++) {
                sizes[i] = split_points_on_sequence[i].size();
            }
        }

        // Broadcast the sizes vector to all ranks
        size_t sizes_size = sizes.size();
        MPI_Bcast(&sizes_size, 1, MPI_UNSIGNED_LONG, 0, MPI_COMM_WORLD);
        sizes.resize(sizes_size);
        MPI_Bcast(sizes.data(), sizes.size(), MPI_UNSIGNED_LONG, 0, MPI_COMM_WORLD);

        // Broadcast the flattened data to all ranks
        size_t total_size = flattened_data.size();
        MPI_Bcast(&total_size, 1, MPI_UNSIGNED_LONG, 0, MPI_COMM_WORLD);  // Envia o tamanho total do vetor
        flattened_data.resize(total_size);
        MPI_Bcast(flattened_data.data(), total_size * sizeof(std::pair<int_t, int_t>), MPI_BYTE, 0, MPI_COMM_WORLD);

        // Unflatten the data received back into the original format
        std::vector<std::vector<std::pair<int_t, int_t>>> received_split_points;
        received_split_points.resize(sizes.size());

        size_t idx = 0;
        for (size_t i = 0; i < sizes.size(); i++) {
            received_split_points[i].resize(sizes[i]);
            for (size_t j = 0; j < sizes[i]; j++) {
                received_split_points[i][j] = flattened_data[idx++];
            }
        }

        split_and_parallel_align(data, name, received_split_points, world_rank, world_size);
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
        output = "(FMAlign2) - total time: " + s.str() + " seconds.";
        print_table_line(output, world_rank);
    }

    // Finalizar MPI
    MPI_Finalize();

    return 0;
}
