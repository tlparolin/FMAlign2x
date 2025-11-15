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

// FMAlign2x - An extended version of FMAlign2 for aligning multiple ultra-long sequences
// Author: Thiago Luiz Parolin
// Contact: thiago.parolin@unesp.br
// July 2025

#ifndef SEQUENCE_SPLIT_ALIGN_H
#define SEQUENCE_SPLIT_ALIGN_H

#include "../ext/SW/ssw.h"
#include "../ext/SW/ssw_cpp.h"
#include "common.h"
#include "utils.h"
#include <algorithm>
#include <format>
#include <future>
#include <iostream>
#include <spoa/spoa.hpp>
#include <sstream>
#include <string>
#include <vector>
#ifdef __linux__
#include <sys/stat.h>
#else
#include <direct.h>
#endif
#include <algorithm>
#include <climits>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <map>
#include <memory>
#include <numeric>
#include <random>
#include <sstream>
#include <string>
#include <string_view>
#include <thread>
#include <tuple>
#include <utility>
#include <vector>

const std::string TMP_FOLDER = "./temp/";

struct ExpandChainParams {
    std::vector<std::string> *data;
    std::vector<std::vector<std::pair<int_t, int_t>>> *chain;
    uint_t chain_index;
    std::vector<std::vector<std::string>>::iterator result_store;
};

struct ParallelAlignParams {
    std::vector<std::string> *data;
    std::vector<std::vector<std::pair<int_t, int_t>>>::iterator parallel_range;
    uint_t task_index;
    std::vector<std::vector<std::string>>::iterator result_store;
    const std::vector<bool> *fallback_needed;
};

struct SpoaTaskParams {
    const std::vector<std::string> *data = nullptr;              // Pointer to input sequences (optional)
    const std::vector<std::pair<int_t, int_t>> *range = nullptr; // Input range (for direct block alignment)
    std::vector<std::string> local_sequences;                    // Used when aligning subdivided sub-blocks
    std::shared_ptr<std::vector<std::string>> result_local;      // Result container for local sub-blocks
    std::vector<std::string> *result_store = nullptr;            // Destination for aligned output
    uint_t seq_num = 0;                                          // Number of sequences in alignment
    uint_t task_index = 0;                                       // Block index
    bool use_batch;                                              // Whether to use batch SPOA alignment
    size_t batch_size;                                           // Batch size for SPOA alignment
};

struct SubBlockInfo {
    size_t seq_id;
    size_t start; // posição no fragmento
    size_t end;   // posição no fragmento
};

void concat_alignment_from_blocks(const std::vector<std::vector<std::string>> &blocks, const std::vector<std::string> &names);

/**
 * @brief Split and parallel align multiple sequences using a vector of chain pairs.
 * This function takes in three parameters: a vector of input sequences (data), a vector of sequence names (name),
 * and a vector of chain pairs (chain) that represent initial pairwise alignments between sequences.
 * It then splits the chain pairs into smaller regions and performs parallel sequence alignment on these regions.
 * Finally, it concatenates the aligned regions and performs sequence-to-profile alignment to generate a final alignment.
 * @param data A vector of input sequences to be aligned
 * @param name A vector of sequence names corresponding to the input sequences
 * @param chain A vector of chain pairs representing initial pairwise alignments between sequences
 * @return void
 */
void split_and_parallel_align(std::vector<std::string> data, std::vector<std::string> name,
                              std::vector<std::vector<std::pair<int_t, int_t>>> chain, ThreadPool &pool);

/**
 * @brief Runs SPOA alignment on a set of sequences (helper for subdivisions)
 * @param seqs Input sequences (may contain empty strings)
 * @return Aligned sequences with gaps
 */
std::vector<std::string> run_spoa_local(const std::vector<std::string> &seqs);

/**
 * @brief Preprocesses alignment blocks between MEMs to reduce load on external aligners.
 * This function attempts to resolve alignment blocks (typically between MEMs) in a fast and memory-efficient
 * way, before falling back to more expensive external aligners such as MAFFT or HAlign. It operates by
 * analyzing each block defined in `parallel_align_range` and choosing one of three strategies:
 * 1. **Exact Match**: If all sequences in the block are identical, simply copies the fragment.
 * 2. **SPOA Alignment**: If average fragment length is small (< 15000 bp), applies SPOA for fast in-memory MSA.
 * 3. **Fallback Flag**: For long or divergent blocks, defers to external aligners by setting a fallback flag.
 * @param data The original vector of sequences (one string per sequence).
 * @param parallel_align_range A vector of alignment ranges (start, length pairs) per sequence, per block.
 *        Each element defines the intervals to be extracted from `data` for one alignment block.
 * @return A pair:
 * - First: A 2D vector of strings containing aligned fragments (SPOA-aligned or copied directly).
 * - Second: A boolean vector where `true` indicates the block must be aligned by an external aligner.
 * @note This function assumes that `spoa_align()` is implemented and available in scope.
 *       It is intended to be used directly after `get_parallel_align_range()` in the FMAlign2 pipeline.
 */
std::vector<std::vector<std::string>>
preprocess_parallel_blocks(const std::vector<std::string> &data,
                           const std::vector<std::vector<std::pair<int_t, int_t>>> &parallel_align_range, ThreadPool &pool);

/**
@brief Expands the chain at the given index for all sequences in the input data.
This function takes a void pointer to input arguments and casts it to the correct struct type.
It then retrieves the required variables, which include the data and chain input parameters, and the chain index.
The query sequence is obtained from the chain for each sequence in the data input, and a reference sequence is obtained
using the neighboring chains. If there is no neighboring chain, the reference sequence is obtained from the start of
the sequence to the end of the previous neighboring chain or to the end of the sequence if there is no previous neighboring chain.
The aligner is then used to align the query sequence and the reference sequence, and the result is stored in an alignment object.
If the chain at the given index for the current sequence is empty, the result is stored in the chain. Otherwise, the query sequence
is already aligned, and the aligned fragment is stored in the aligned_fragment vector.
Finally, the function stores the aligned fragments in the result_store vector.
@param arg A void pointer to input arguments.
@return NULL
*/
void *expand_chain(void *arg);

/**
 * @brief Store the Smith-Waterman alignment results in a vector of aligned sequences.
 * This function takes the alignment results generated by the StripedSmithWaterman algorithm and
 * stores the aligned sequences in a vector of strings. The function also returns the alignment
 * start and length as a pair of integers. If the alignment failed, the function returns (-1,-1).
 * @param alignment The alignment results generated by StripedSmithWaterman algorithm.
 * @param ref The reference sequence.
 * @param query The query sequence.
 * @param res_store A vector of aligned sequences.
 * @param seq_index The index of the query sequence in the vector of aligned sequences.
 * @return A pair of integers representing the alignment start and length.
 */
std::pair<int_t, int_t> store_sw_alignment(StripedSmithWaterman::Alignment alignment, std::string &ref, std::string &query,
                                           std::vector<std::string> &res_store, uint_t seq_index);

/**
 * @brief Get the range of each sequence in parallel alignment
 * @param data The vector of sequences to be aligned
 * @param chain The vector of chains representing the alignment
 * @return The vector of ranges for each sequence in the alignment
 */
std::vector<std::vector<std::pair<int_t, int_t>>> get_parallel_align_range(const std::vector<std::string> &data,
                                                                           const std::vector<std::vector<std::pair<int_t, int_t>>> &chain);

/**
 * @brief Concatenate two sets of sequence data (chain and parallel) into a single set of concatenated data.
 * @param chain_string A vector of vectors containing the chain sequence data.
 * @param parallel_string A vector of vectors containing the parallel sequence data.
 * @return std::vector<std::vectorstd::string> A vector of vectors containing the concatenated sequence data.
 */
std::vector<std::vector<std::string>> concat_chain_and_parallel(const std::vector<std::vector<std::string>> &chain_string,
                                                                const std::vector<std::vector<std::string>> &parallel_string);

#endif