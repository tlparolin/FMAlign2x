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
// Nov 2025

#include "sequence_split_align.h"

/**
 * @brief Generates a random string of the specified length.
 * This function generates a random string of the specified length. The generated string
 * consists of lowercase English letters ('a' to 'z'). It uses a thread-local random number generator.
 * @param length The length of the generated string.
 * @return A random string of the specified length, or an empty string if an error occurs.
 */
std::string generateRandomString(int length) {
    // Conjunto de caracteres permitido
    constexpr std::string_view charset = "abcdefghijklmnopqrstuvwxyz";
    // Garante eficiência, pois evita reallocs
    std::string result;
    result.reserve(length);

    // Geradores thread-safe
    static thread_local std::mt19937 gen{std::random_device{}()};
    static thread_local std::uniform_int_distribution<std::size_t> dist(0, charset.size() - 1);

    // Algoritmo moderno, funcional e rápido
    std::ranges::generate_n(std::back_inserter(result), length, [&] { return charset[dist(gen)]; });

    return result;
}

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
std::string random_file_end;

void split_and_parallel_align(std::vector<std::string> data, std::vector<std::string> name,
                              std::vector<std::vector<std::pair<int_t, int_t>>> chain, ThreadPool &pool) {
    if (global_args.verbose) {
        std::cout << "#                Parallel Aligning...                       #" << std::endl;
        print_table_divider();
    }

    random_file_end = generateRandomString(10);
    std::string output = "";
    Timer timer;
    uint_t chain_num = chain[0].size();
    uint_t seq_num = data.size();
    std::vector<std::vector<std::string>> chain_string(chain_num);

    // === EXPAND CHAINS (Smith-Waterman local) ===
    std::vector<ExpandChainParams> params(chain_num);
    for (uint_t i = 0; i < chain_num; i++) {
        params[i].data = &data;
        params[i].chain = &chain;
        params[i].chain_index = i;
        params[i].result_store = chain_string.begin() + i;
    }

    if (global_args.min_seq_coverage == 1) {
        for (uint_t i = 0; i < chain_num; i++) {
            pool.add_task([i, &params]() { expand_chain(&params[i]); });
        }
        pool.wait_for_tasks();
    } else {
        for (uint_t i = 0; i < chain_num; i++) {
            expand_chain(&params[i]);
        }
    }

    params.clear();

    if (global_args.verbose) {
        print_table_line(std::format("SW expand time: {:.2f} seconds.", timer.elapsed_time()));
    }

    timer.reset();

    // === PREPARE PARALLEL ALIGNMENT RANGES ===
    std::vector<std::vector<std::pair<int_t, int_t>>> parallel_align_range = get_parallel_align_range(data, chain);

    // === SPOA-BASED ALIGNMENT (IN-MEMORY ONLY) ===
    std::vector<std::vector<std::string>> spoa_parallel_string;

    spoa_parallel_string = preprocess_parallel_blocks(data, parallel_align_range, pool);

    if (global_args.verbose) {
        print_table_line(std::format("Parallel Align Time: {:.2f} seconds.", timer.elapsed_time()));
    }

    timer.reset();
    // fallback to empty structure if not extended
    spoa_parallel_string.resize(parallel_align_range.size(), std::vector<std::string>(seq_num, ""));

    // === CONCATENATION ===
    std::vector<std::vector<std::pair<int_t, int_t>>> concat_range = concat_chain_and_parallel_range(chain, parallel_align_range);

    // Normalize blocks that are empty
    for (auto &block : spoa_parallel_string) {
        if (block.size() != seq_num) {
            size_t L = 0;
            for (const auto &s : block)
                L = std::max(L, s.size());
            block.resize(seq_num);
            std::fill(block.begin() + (block.size() - (seq_num - block.size())), block.end(), std::string(L, '-'));
        }
    }

    std::vector<std::vector<std::string>> concat_string = concat_chain_and_parallel(chain_string, spoa_parallel_string);
    std::vector<uint_t> fragment_len = get_first_nonzero_lengths(concat_string);

    seq2profile(concat_string, data, concat_range, fragment_len);

    concat_alignment(concat_string, name);

    if (global_args.verbose) {
        print_table_line(std::format("Seq-profile time: {:.2f} seconds.", timer.elapsed_time()));
        print_table_divider();
    }
}

/**
 * @brief Executes a SPOA (Simd Partial Order Alignment) task on a set of DNA sequence fragments.
 * This function is designed to be run as a thread. It extracts subsequences from the input data
 * based on the provided ranges, constructs a vector of fragments, and performs multiple sequence
 * alignment using the SPOA algorithm. The result is stored in the location pointed to by
 * `params->result_store`.
 * @param arg A pointer to a `SpoaTaskParams` structure containing:
 *  - the number of sequences (`seq_num`),
 *  - the input data (`data`),
 *  - the ranges to extract from each sequence (`range`),
 *  - and a pointer to where the alignment result should be stored (`result_store`).
 * @return Always returns `nullptr`. The result of the alignment is stored via the pointer in `params`.
 * @note If a sequence has an invalid range (start == -1 or length <= 0), an empty string is used.
 */
void *spoa_task(void *arg) {
    // Cast the argument to the expected structure
    SpoaTaskParams *params = static_cast<SpoaTaskParams *>(arg);

    // Prepare a vector to hold extracted fragments for each sequence
    std::vector<std::string> fragments(params->seq_num);

    // Loop over each sequence to extract the corresponding fragment based on start and length
    for (uint_t s = 0; s < params->seq_num; ++s) {
        auto [start, len] = (*params->range)[s];
        if (start != -1 && len > 0) {
            // Extract a valid fragment from sequence
            fragments[s] = (*params->data)[s].substr(start, len);
        } // } else {
        //     // For invalid range, set fragment as empty string
        //     fragments[s] = "";
        // }
    }

    // Run the SPOA alignment on the extracted fragments
    // and store the result at the designated location
    *(params->result_store) = spoa_align(fragments);

    // Return nullptr as this function is intended for threading interface compatibility
    return nullptr;
}

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
                           const std::vector<std::vector<std::pair<int_t, int_t>>> &parallel_align_range, ThreadPool &pool) {
    const size_t MAX_BLOCK_SIZE = 15000; // Maximum block size for SPOA (in bases)
    const size_t OVERLAP_SIZE = 200;     // Overlap between consecutive sub-blocks

    uint_t parallel_num = parallel_align_range.size();
    uint_t seq_num = data.size();

    std::vector<std::vector<std::string>> fast_parallel_string(parallel_num, std::vector<std::string>(seq_num));

    uint_t count_exact = 0;
    uint_t count_spoa = 0;
    uint_t count_split = 0;

    // Structure to hold subdivision info
    struct SubBlockInfo {
        uint_t block_index;
        std::vector<std::vector<std::string>> sub_results; // One result per sub-block
        std::vector<size_t> sub_block_starts;              // Start position of each sub-block
        std::vector<size_t> sub_block_ends;                // End position of each sub-block
    };

    std::vector<SubBlockInfo> subdivided_blocks;
    std::vector<SpoaTaskParams> spoa_params;

    // ============================================================
    // PASS 1: Classify blocks and prepare tasks
    // ============================================================
    for (uint_t i = 0; i < parallel_num; ++i) {
        const auto &range = parallel_align_range[i];
        std::vector<std::string> fragments(seq_num);
        bool all_equal = true;
        size_t total_len = 0;
        size_t max_len = 0;

        // Extract fragments
        for (uint_t s = 0; s < seq_num; ++s) {
            auto [start, len] = range[s];
            if (start != -1 && len > 0) {
                fragments[s] = data[s].substr(start, len);
                total_len += len;
                max_len = std::max(max_len, fragments[s].size());
                if (s > 0 && fragments[s] != fragments[0])
                    all_equal = false;
            } else {
                fragments[s] = "";
            }
        }

        size_t avg_len = (seq_num > 0) ? (total_len / seq_num) : 0;

        // --- Case 1: All fragments identical → direct copy ---
        if (all_equal) {
            std::fill(fast_parallel_string[i].begin(), fast_parallel_string[i].end(), fragments[0]);
            ++count_exact;
            continue;
        }

        // --- Case 2: Small block → direct SPOA ---
        if (avg_len <= MAX_BLOCK_SIZE) {
            SpoaTaskParams params;
            params.data = &data;
            params.range = &parallel_align_range[i];
            params.task_index = i;
            params.seq_num = seq_num;
            params.result_store = &fast_parallel_string[i];
            params.result_local = nullptr; // Not used for direct SPOA
            spoa_params.push_back(params);
            ++count_spoa;
            continue;
        }

        // --- Case 3: Large block → subdivide ---
        ++count_split;

        SubBlockInfo sub_info;
        sub_info.block_index = i;

        size_t pos = 0;

        // Generate overlapping sub-blocks
        while (pos < max_len) {
            size_t end = std::min(pos + MAX_BLOCK_SIZE, max_len);

            // Add overlap to the end (except for last block)
            if (end < max_len) {
                end = std::min(end + OVERLAP_SIZE, max_len);
            }

            sub_info.sub_block_starts.push_back(pos);
            sub_info.sub_block_ends.push_back(end);

            // Extract sub-fragments
            std::vector<std::string> sub_fragments(seq_num);
            for (uint_t s = 0; s < seq_num; ++s) {
                const std::string &seq = fragments[s];
                if (seq.empty()) {
                    sub_fragments[s] = "";
                } else {
                    size_t real_start = std::min(pos, seq.size());
                    size_t real_end = std::min(end, seq.size());
                    if (real_start < real_end) {
                        sub_fragments[s] = seq.substr(real_start, real_end - real_start);
                    } else {
                        sub_fragments[s] = "";
                    }
                }
            }

            // Create SPOA task for this sub-block
            SpoaTaskParams sub_params;
            sub_params.task_index = i;
            sub_params.seq_num = seq_num;
            sub_params.data = nullptr;
            sub_params.range = nullptr;
            sub_params.result_store = nullptr;
            sub_params.local_sequences = sub_fragments;
            sub_params.result_local = std::make_shared<std::vector<std::string>>();
            spoa_params.push_back(sub_params);
            ++count_spoa;

            // Move to next sub-block (with overlap)
            pos += (MAX_BLOCK_SIZE - OVERLAP_SIZE);
            if (pos >= max_len)
                break;
        }

        subdivided_blocks.push_back(std::move(sub_info));
    }

    // ============================================================
    // PASS 2: Execute all SPOA tasks in parallel
    // ============================================================
    if (!spoa_params.empty()) {

        for (auto &params : spoa_params) {
            pool.add_task([&params]() {
                if (params.data != nullptr) {
                    // Direct SPOA task (small block)
                    spoa_task(&params);
                } else {
                    // Subdivision SPOA task
                    std::vector<std::string> aligned = run_spoa_local(params.local_sequences);
                    *(params.result_local) = aligned;
                }
            });
        }

        pool.wait_for_tasks();
    }

    // ============================================================
    // PASS 3: Merge subdivided blocks with overlap trimming
    // ============================================================
    static const std::vector<std::vector<std::string>> empty_results;

    // Populate a map for quick access to sub-block results
    std::unordered_map<uint_t, std::vector<std::vector<std::string>>> block_results_map;
    for (const auto &p : spoa_params) {
        if (p.result_local != nullptr) {
            block_results_map[p.task_index].push_back(*(p.result_local));
        }
    }

    // Process each subdivided block
    for (auto &sub_info : subdivided_blocks) {
        uint_t block_idx = sub_info.block_index;
        std::vector<std::string> merged(seq_num, "");

        // Collect all sub-block results
        auto it = block_results_map.find(block_idx);
        const auto &sub_results = (it != block_results_map.end()) ? it->second : empty_results;

        if (sub_results.empty()) {
            // Safety: if no results, fill with gaps
            size_t expected_len = 0;
            for (uint_t s = 0; s < seq_num; ++s) {
                auto [start, len] = parallel_align_range[block_idx][s];
                expected_len = std::max(expected_len, static_cast<size_t>(len));
            }
            fast_parallel_string[block_idx] = std::vector<std::string>(seq_num, std::string(expected_len, '-'));
            continue;
        }

        // Merge sub-blocks
        std::vector<std::string> overlap_prev(seq_num);
        std::vector<std::string> overlap_curr(seq_num);
        std::vector<std::string> bridge_fragments(seq_num);

        for (size_t sub_idx = 0; sub_idx < sub_results.size(); ++sub_idx) {
            const auto &sub_result = sub_results[sub_idx];

            if (sub_idx == 0) {
                // First sub-block: copy everything
                for (uint_t s = 0; s < seq_num; ++s) {
                    merged[s] = sub_result[s];
                }
            } else {
                // Subsequent sub-blocks: trim overlap and concatenate
                const size_t overlap_in_coords = OVERLAP_SIZE;

                // Extract the top of the previous block (merged) and the beginning of the current block (sub_result)

                for (size_t s = 0; s < seq_num; ++s) {
                    const std::string &prev_seq = merged[s];
                    const std::string &curr_seq = sub_result[s];

                    size_t len_prev = std::min(overlap_in_coords, prev_seq.size());
                    size_t len_curr = std::min(overlap_in_coords, curr_seq.size());

                    overlap_prev[s] = prev_seq.substr(prev_seq.size() - len_prev);
                    overlap_curr[s] = curr_seq.substr(0, len_curr);
                }

                // Creates intermediate sequence by joining edges
                for (size_t s = 0; s < seq_num; ++s)
                    bridge_fragments[s] = overlap_prev[s] + overlap_curr[s];

                // Align the bridge with SPOA
                std::vector<std::string> bridge_aln = run_spoa_local(bridge_fragments);

                // Choose cut point in the middle of the bridge alignment
                size_t bridge_cut = bridge_aln[0].size() / 2;

                // Apply SPOA merge of refined edges
                for (size_t s = 0; s < seq_num; ++s) {
                    size_t trim_len = std::min(overlap_in_coords, merged[s].size());
                    merged[s].erase(merged[s].size() - trim_len);  // remove fim antigo
                    merged[s] += bridge_aln[s].substr(bridge_cut); // anexa parte refinada
                }
            }
        }

        fast_parallel_string[block_idx] = merged;
    }

    // ============================================================
    // Verbose output
    // ============================================================
    if (global_args.verbose) {
        print_table_line("Blocks copied directly: " + std::to_string(count_exact));
        print_table_line("Blocks aligned by SPOA: " + std::to_string(count_spoa));
        print_table_line("Large blocks subdivided: " + std::to_string(count_split));
    }

    return fast_parallel_string;
}

/**
 * @brief Runs SPOA alignment on a set of sequences (helper for subdivisions)
 * @param seqs Input sequences (may contain empty strings)
 * @return Aligned sequences with gaps
 */
std::vector<std::string> run_spoa_local(const std::vector<std::string> &seqs) {
    if (seqs.empty())
        return {};

    // Use SPOA parameters (match=5, mismatch=-4, gap=-8)
    auto aligner = spoa::AlignmentEngine::Create(spoa::AlignmentType::kNW, 5, -4, -8);
    spoa::Graph graph;

    std::vector<size_t> non_empty_indices;
    non_empty_indices.reserve(seqs.size());

    // Add only non-empty sequences
    for (size_t i = 0; i < seqs.size(); ++i) {
        if (seqs[i].empty())
            continue;

        auto alignment = aligner->Align(seqs[i], graph);
        graph.AddAlignment(alignment, seqs[i]);
        non_empty_indices.push_back(i);
    }

    // If all empty, return empty strings
    if (non_empty_indices.empty()) {
        return std::vector<std::string>(seqs.size(), "");
    }

    // Generate MSA
    auto msa_compact = graph.GenerateMultipleSequenceAlignment();
    size_t aln_len = msa_compact.empty() ? 0 : msa_compact[0].size();

    // Rebuild to original size
    std::vector<std::string> msa(seqs.size(), std::string(aln_len, '-'));
    for (size_t k = 0; k < non_empty_indices.size(); ++k) {
        msa[non_empty_indices[k]] = std::move(msa_compact[k]);
    }

    return msa;
}

/**
 * @brief Performs multiple sequence alignment using SPOA (Partial Order Alignment).
 * This function leverages the SPOA library to construct a partial order graph and
 * progressively align input sequences to it. It uses the global alignment model
 * (Needleman-Wunsch) with linear gap penalties, optimized for short to medium-length
 * sequence blocks (e.g., between MEMs). Empty sequences are ignored during the alignment process.
 * The final output is a multiple sequence alignment with gaps introduced as necessary to maintain consistency.
 * Scoring parameters used:
 * - Match: +4
 * - Mismatch: -10
 * - Gap (linear): -8
 * @param sequences A vector of input sequences to be aligned. Each sequence is a std::string.
 * @return A vector of aligned sequences (same size as input), each padded with '-' where needed.
 * @note This function is designed for fast in-memory alignment and is suitable as a lightweight
 * alternative to full external aligners (e.g., MAFFT, HAlign) in internal blocks of ultralong sequences.
 * @see https://github.com/rvaser/spoa for more details on the SPOA library.
 */
std::vector<std::string> spoa_align(const std::vector<std::string> &sequences) {
    if (sequences.empty())
        return {};

    // Create the SPOA alignment engine with linear gap penalties
    // Note: Each thread creates its own engine instance for thread safety
    auto alignment_engine = spoa::AlignmentEngine::Create(spoa::AlignmentType::kNW, 5, -4, -8);

    spoa::Graph graph{};

    // maps indices of non-empty sequences
    std::vector<size_t> map_idx;
    map_idx.reserve(sequences.size());
    for (size_t i = 0; i < sequences.size(); ++i) {
        const auto &seq = sequences[i];
        if (seq.empty())
            continue;
        auto aln = alignment_engine->Align(seq, graph);
        graph.AddAlignment(aln, seq);
        map_idx.push_back(i);
    }

    // If all are empty: return N empty strings
    if (map_idx.empty()) {
        return std::vector<std::string>(sequences.size(), std::string{});
    }

    // Generates MSA only from non-empty ones
    auto msa_compact = graph.GenerateMultipleSequenceAlignment();
    const size_t aln_len = msa_compact.empty() ? 0 : msa_compact.front().size();

    // Rebuilds to original size: empty spaces become just gaps
    std::vector<std::string> msa(sequences.size(), std::string(aln_len, '-'));
    for (size_t k = 0; k < msa_compact.size(); ++k) {
        msa[map_idx[k]] = std::move(msa_compact[k]);
    }
    return msa;
}

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
void *expand_chain(void *arg) {
    auto *ptr = static_cast<ExpandChainParams *>(arg);
    const auto &data = *(ptr->data); // No copy, only reference
    auto &chain = *(ptr->chain);     // Reference, allows mutation
    const auto chain_index = ptr->chain_index;
    const size_t seq_num = data.size();
    const size_t chain_num = chain[0].size();

    StripedSmithWaterman::Aligner aligner;
    StripedSmithWaterman::Filter filter;
    StripedSmithWaterman::Alignment alignment;

    // Find first valid query for the current chain index
    std::string query;
    int_t query_length = 0;
    bool found_query = false;
    for (size_t i = 0; i < seq_num; ++i) {
        if (chain[i][chain_index].first != -1) {
            query_length = chain[i][chain_index].second;
            query = data[i].substr(chain[i][chain_index].first, query_length);
            found_query = true;
            break;
        }
    }

    if (!found_query) {
        // No valid query found; result is empty, safe fallback
        *(ptr->result_store) = std::vector<std::string>(seq_num, "");
        return nullptr;
    }

    // Prepare vector with fixed size for results
    std::vector<std::string> aligned_fragment(seq_num);

    for (size_t i = 0; i < seq_num; ++i) {
        int_t begin_pos = chain[i][chain_index].first;
        if (begin_pos == -1) {
            size_t tmp_index = chain_index;
            int_t maskLen = std::max(query_length / 2, 15);

            // Find reference boundaries for alignment
            size_t ref_begin_pos = 0;
            if (tmp_index > 0 && chain[i][tmp_index - 1].first != -1) {
                ref_begin_pos = chain[i][tmp_index - 1].first + chain[i][tmp_index - 1].second;
            }
            size_t ref_end_pos = data[i].length() - 1;
            if (tmp_index < chain_num - 1 && chain[i][tmp_index + 1].first != -1) {
                ref_end_pos = chain[i][tmp_index + 1].first;
            }
            std::string ref = data[i].substr(ref_begin_pos, ref_end_pos - ref_begin_pos);

            // Perform alignment
            aligner.Align(query.c_str(), ref.c_str(), ref.size(), filter, &alignment, maskLen);
            auto p = store_sw_alignment(alignment, ref, query, aligned_fragment, i);
            if (p.first != -1) {
                p.first += ref_begin_pos;
                chain[i][chain_index] = p; // update chain with alignment info
            }
        } else {
            // Already aligned: copy query
            aligned_fragment[i] = query;
        }
    }

    *(ptr->result_store) = std::move(aligned_fragment);
    return nullptr;
}

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
                                           std::vector<std::string> &res_store, uint_t seq_index) {
    // Get a constant reference to the cigar vector from the alignment to avoid copying
    std::vector<unsigned int> const &cigar = alignment.cigar;
    int_t ref_begin = alignment.ref_begin;
    int_t ref_end = alignment.ref_end;
    uint_t query_begin = 0;
    // If the alignment failed (ref_begin < 0), clear the result and return failure pair (-1, -1)
    if (ref_begin < 0) {
        res_store[seq_index].clear();
        return {-1, -1};
    }

    // Variables to count soft clipped bases and total query length
    int_t S_count = 0;
    int_t total_length = alignment.query_end;
    int_t new_ref_begin = ref_begin;
    uint_t new_ref_end = ref_end;

    // Lambdas to extract length and operation from cigar integer encoding for convenience
    auto cigar_len = [](auto c) { return cigar_int_to_len(c); };
    auto cigar_op = [](auto c) { return cigar_int_to_op(c); };

    // Adjust new_ref_begin and S_count if the first cigar operation is soft clip (S)
    if (cigar_op(cigar.front()) == 'S') {
        auto len = cigar_len(cigar.front());
        S_count += len;
        new_ref_begin = std::max(int_t{0}, new_ref_begin - static_cast<int_t>(len));
    }
    // Adjust new_ref_end and S_count if the last cigar operation is soft clip (S)
    if (cigar_op(cigar.back()) == 'S') {
        auto len = cigar_len(cigar.back());
        S_count += len;
        total_length += len;
        new_ref_end = std::min(uint_t(ref.length() - 1), new_ref_end + len);
    }

    // If the amount of soft clipping is too large relative to the total length, clear result and return failure
    if (S_count > static_cast<int_t>(std::ceil(0.8 * total_length))) {
        res_store[seq_index].clear();
        return {-1, -1};
    }

    // String to accumulate the aligned reference sequence result
    std::string aligned_result;
    // Pointers tracking positions in reference and query sequences
    uint_t p_ref = 0, p_query = 0;

    // Iterate over each cigar operation, handling accordingly
    for (auto const &c : cigar) {
        char op = cigar_op(c);
        int_t len = cigar_len(c);
        switch (op) {
        case 'S': {
            // Handle soft clipping differently for the start and end of the alignment
            if (&c == &cigar.front()) {
                int_t tmp_len = len;
                // If clipping length exceeds beginning of alignment, pad with gaps
                if (ref_begin <= len) {
                    while (tmp_len > ref_begin) {
                        aligned_result += '-';
                        --tmp_len;
                    }
                }
                // Append clipped reference bases after gaps
                for (int_t j = ref_begin - tmp_len; j < ref_begin; ++j) {
                    aligned_result += ref[j];
                }
            } else {
                int_t tmp_len = len;
                // Append clipped reference bases immediately after alignment end
                for (uint_t j = ref_end + 1; j < ref.length() && tmp_len > 0; ++j) {
                    aligned_result += ref[j];
                    --tmp_len;
                }
                // Pad with gaps if clipped length exceeds available reference
                while (tmp_len > 0) {
                    aligned_result += '-';
                    --tmp_len;
                }
            }
            p_query += len; // Advance query pointer by soft clip length
            break;
        }
        case 'M':
        case 'X':
        case '=': {
            // For match, mismatch and equal operations, append reference substring of length len
            aligned_result.append(ref.begin() + ref_begin + p_ref, ref.begin() + ref_begin + p_ref + len);
            p_ref += len;   // Advance reference pointer
            p_query += len; // Advance query pointer
            break;
        }
        case 'I': {
            // For insertion operations, append gap characters to alignment
            aligned_result.append(len, '-');
            p_query += len; // Advance query pointer
            break;
        }
        case 'D': {
            // For deletions, append deleted reference bases
            std::string gaps(len, '-'); // Create gap string to insert into query sequences
            aligned_result.append(ref.begin() + ref_begin + p_ref, ref.begin() + ref_begin + p_ref + len);
            p_ref += len; // Advance reference pointer

            // Insert gaps into query and all previous result sequences at correct positions
            query.insert(query_begin + p_query, gaps);
            for (uint_t j = 0; j < seq_index; ++j) {
                if (!res_store[j].empty()) {
                    res_store[j].insert(query_begin + p_query, gaps);
                }
            }
            p_query += len; // Advance query pointer for deletion length
            break;
        }
        }
    }

    // Store the aligned result string in the results vector at the appropriate index
    res_store[seq_index] = std::move(aligned_result);

    // Return the adjusted reference begin position and length of the alignment on the reference
    return {new_ref_begin, new_ref_end - new_ref_begin + 1};
}

/**
 * @brief Get the range of each sequence in parallel alignment
 * @param data The vector of sequences to be aligned
 * @param chain The vector of chains representing the alignment
 * @return The vector of ranges for each sequence in the alignment
 */
std::vector<std::vector<std::pair<int_t, int_t>>> get_parallel_align_range(const std::vector<std::string> &data,
                                                                           const std::vector<std::vector<std::pair<int_t, int_t>>> &chain) {
    uint_t seq_num = data.size();
    // Get maximum size for sequences
    uint_t chain_num =
        std::max_element(chain.begin(), chain.end(), [](const auto &a, const auto &b) { return a.size() < b.size(); })->size();

    // Pre-allocate memory to avoid resizing during push_back
    std::vector<std::vector<std::pair<int_t, int_t>>> parallel_align_range(seq_num);
    parallel_align_range.reserve(seq_num); // Reserve space to avoid reallocations

    // Iterate through each sequence
    for (uint_t i = 0; i < seq_num; i++) {
        int_t last_pos = 0;
        std::vector<std::pair<int_t, int_t>> &tmp_range = parallel_align_range[i];
        tmp_range.reserve(chain_num + 1); // Reserve enough space for the ranges

        // Iterate through each chain for the current sequence
        for (uint_t j = 0; j < chain_num; j++) {
            int_t begin_pos = chain[i][j].first;

            if (begin_pos == -1) {
                tmp_range.push_back({-1, -1});
                last_pos = -1;
            } else {
                if (last_pos == -1) {
                    tmp_range.push_back({-1, -1});
                } else {
                    tmp_range.push_back({last_pos, begin_pos - last_pos});
                }
                last_pos = begin_pos + chain[i][j].second;
            }
        }

        // Handle the end of the sequence after the last chain
        if (last_pos == -1) {
            tmp_range.push_back({-1, -1});
        } else {
            tmp_range.push_back({last_pos, data[i].length() - last_pos});
        }
    }

    // Allocate transpose result vector
    std::vector<std::vector<std::pair<int_t, int_t>>> transpose_res;
    transpose_res.reserve(chain_num + 1);

    // Transpose parallel_align_range into transpose_res
    for (uint_t j = 0; j <= chain_num; j++) {
        std::vector<std::pair<int_t, int_t>> transposed_row(seq_num);
        for (uint_t i = 0; i < seq_num; i++) {
            transposed_row[i] = parallel_align_range[i][j];
        }
        transpose_res.push_back(std::move(transposed_row)); // Move to avoid copies
    }

    return transpose_res;
}

/**
 * @brief Concatenate multiple sequence alignments into a single alignment and write the result to an output file.
 * @param concat_string A 2D vector of strings containing the aligned sequences to concatenate.
 * @param name A vector of strings containing the names of the sequences.
 */
void concat_alignment(const std::vector<std::vector<std::string>> &concat_string, const std::vector<std::string> &name) {
    std::string output_path = global_args.output_path;
    std::vector<std::string> concated_data(name.size());

    // Estimate final size and reserve buffer for each sequence
    for (size_t i = 0; i < name.size(); ++i) {
        size_t total_len = 0;
        for (size_t j = 0; j < concat_string.size(); ++j) {
            total_len += concat_string[j][i].size();
        }
        concated_data[i].reserve(total_len);
    }

    // Concatenate fragments efficiently
    for (size_t i = 0; i < name.size(); ++i) {
        for (size_t j = 0; j < concat_string.size(); ++j) {
            concated_data[i].append(concat_string[j][i]);
        }
    }

    // Write output in FASTA format
    std::ofstream output_file(output_path);
    if (!output_file) {
        throw std::runtime_error("Error opening output file " + output_path);
    }

    for (size_t i = 0; i < concated_data.size(); ++i) {
        output_file << ">" << name[i] << "\n" << concated_data[i] << "\n";
    }
    // output_file closes automatically by RAII
}

bool cmp(const std::pair<uint_t, uint_t> &a, const std::pair<uint_t, uint_t> &b) { return a.second < b.second; }

/**
 * @brief Convert sequence fragments into profile by aligning missing fragments with existing ones.
 * @param concat_string A reference to a vector of vectors of strings representing concatenated sequence fragments.
 * @param data A reference to a vector of strings representing the sequence names.
 * @param concat_range A reference to a vector of vectors of pairs of integers representing the start and end positions of the sequence
 * fragments.
 * @param fragment_len A reference to a vector of unsigned integers representing the lengths of the sequence fragments.
 * @return None.
 */
void seq2profile(std::vector<std::vector<std::string>> &concat_string, std::vector<std::string> &data,
                 std::vector<std::vector<std::pair<int_t, int_t>>> &concat_range, std::vector<uint_t> &fragment_len) {
    // Count the number of missing fragments for each sequence and store them in a vector of pairs.
    uint_t seq_num = data.size();
    uint_t fragment_num = fragment_len.size();

    std::vector<std::pair<uint_t, uint_t>> missing_fragment_count(seq_num);
    for (uint_t i = 0; i < seq_num; i++) {
        uint_t count = 0;
        for (uint_t j = 0; j < fragment_num; j++) {
            if (concat_range[j][i].first == -1) {
                count += fragment_len[j];
            }
        }
        missing_fragment_count[i] = std::make_pair(i, count);
    }

    // Sort the vector of pairs by the number of missing fragments in ascending order.
    std::sort(missing_fragment_count.begin(), missing_fragment_count.end(), cmp);
#if DEBUG
    for (int_t i = 0; i < missing_fragment_count.size(); i++) {
        std::cout << missing_fragment_count[i].first << " " << missing_fragment_count[i].second << std::endl;
    }
#endif //

    // For each sequence with missing fragments, align the missing fragments with existing fragments.
    for (uint_t i = 0; i < seq_num; i++) {
        uint_t seq_index = missing_fragment_count[i].first;
        uint_t missing_count = missing_fragment_count[i].second;
        if (missing_count <= 0) {
            continue;
        }
        // std::cout << seq_index << std::endl;
        bool if_start = false;
        std::vector<std::vector<std::string>>::iterator left_it = concat_string.begin();
        std::vector<std::vector<std::string>>::iterator right_it = concat_string.begin();
        std::vector<std::vector<std::string>>::iterator cur_it = concat_string.begin();
        for (; cur_it != concat_string.end(); cur_it++) {
            std::vector<std::string> cur_vec = *cur_it;
            // std::cout << cur_vec[seq_index].length() << " " << fragment_len[cur_it - concat_string.begin()] << " " << concat_range[cur_it
            // - concat_string.begin()][seq_index].first << std::endl; if (cur_vec[seq_index].length() != fragment_len[cur_it -
            // concat_string.begin()]) {
            if (concat_range[cur_it - concat_string.begin()][seq_index].first == -1) {
                if (!if_start) {
                    left_it = cur_it;
                    if_start = true;
                }
            } else if (if_start) {
                right_it = cur_it - 1;
                cur_it = seq2profile_align(seq_index, left_it - concat_string.begin(), right_it - concat_string.begin(), concat_string,
                                           data, concat_range, fragment_len);
                if_start = false;
            }
            if (if_start && cur_it + 1 == concat_string.end()) {
                right_it = cur_it;
                seq2profile_align(seq_index, left_it - concat_string.begin(), right_it - concat_string.begin(), concat_string, data,
                                  concat_range, fragment_len);
                break;
            }
        }
    }
}

/**
 * @brief: Aligns a sequence and a profile using a third-party tool called profile_two_align and returns the iterator pointing to the next
 * position in the 2D vector of strings. The function first checks for gaps between the fragments and adds them as necessary. Then it writes
 * the sequence to align and a profile file, passes them to profile_two_align, and reads the results. Finally, it updates the concatenated
 * string, range and fragment length information accordingly.
 * @param seq_index: Index of the sequence to align.
 * @param left_index: Index of the left-most fragment.
 * @param right_index: Index of the right-most fragment.
 * @param concat_string: 2D vector of strings containing the concatenated fragments.
 * @param data: Vector of strings containing the sequences.
 * @param concat_range: 2D vector of pairs of integers representing the range of each fragment in each sequence.
 * @param fragment_len: Vector of unsigned integers representing the length of each fragment.
 * @return std::vector<std::vectorstd::string>::iterator: Iterator pointing to the next position in the 2D vector of strings.
 */
std::vector<std::vector<std::string>>::iterator seq2profile_align(uint_t seq_index, uint_t left_index, uint_t right_index,
                                                                  std::vector<std::vector<std::string>> &concat_string,
                                                                  std::vector<std::string> &data,
                                                                  std::vector<std::vector<std::pair<int_t, int_t>>> &concat_range,
                                                                  std::vector<uint_t> &fragment_len) {
    int_t seq_begin = 0;
    // determine the start position of the sequence, if there is a sequence on the left side, take the end of the last one.
    if (left_index >= 1) {
        seq_begin = concat_range[left_index - 1][seq_index].first + concat_range[left_index - 1][seq_index].second;
    }
    int_t seq_end = data[seq_index].length();
    // determine the end position of the sequence, if there is a sequence on the right side, take the start position of the next one.
    if (right_index + 1 < concat_range.size()) {
        seq_end = concat_range[right_index + 1][seq_index].first;
    }
    // create a file to store the sequence content.

    std::string seq_content = data[seq_index].substr(seq_begin, seq_end - seq_begin);
#if DEBUG
    std::cout << seq_content << std::endl;
    std::cout << concat_range.size() << " " << concat_range[0].size() << std::endl;
    std::cout << left_index << " " << right_index << " " << seq_index << std::endl;
    std::cout << seq_end << " " << seq_begin << " " << data[seq_index].length() << std::endl;
#endif // DEBUG

    if (seq_content.length() == 0) {
        for (uint_t i = left_index; i <= right_index; i++) {
            std::string tmp_gaps = "";
            for (uint_t j = 0; j < fragment_len[i]; j++) {
                tmp_gaps += "-";
            }
            concat_string[i][seq_index] = tmp_gaps;
        }

        return concat_string.begin() + right_index + 1;
    }
    std::string seq_file_name = TMP_FOLDER + "seq-" + std::to_string(seq_index) + "_" + random_file_end + ".fasta";
    std::ofstream file;
    file.open(seq_file_name);
    if (!file.is_open()) {
        std::cerr << seq_file_name << " fail to open!" << std::endl;
        exit(1);
    }
    // write the sequence content to the file in fasta format.
    std::stringstream sstream;
    sstream << ">SEQUENCE" << std::to_string(seq_index) << "\n";
    sstream << seq_content << "\n";
    file << sstream.str();
    file.close();

    // create a file to store the profile of the sequences.
    sstream.str("");
    std::string profile_file_name = TMP_FOLDER + "profile-" + std::to_string(seq_index) + "_" + random_file_end + ".fasta";

    file.open(profile_file_name);

    if (!file.is_open()) {
        std::cerr << seq_file_name << " fail to open!" << std::endl;
        exit(1);
    }

    std::vector<uint_t> selected_profile_seq_index;
    std::vector<uint_t> missing_profile_seq_index;
    for (uint_t i = 0; i < data.size(); i++) {
        if (i == seq_index) {
            continue;
        }
        bool cur_fragment_is_missing = false;
        sstream << ">SEQUENCE" << std::to_string(i) << "\n";
        for (uint_t j = left_index; j <= right_index; j++) {
            if (concat_string[j][i].length() == fragment_len[j]) {
                sstream << concat_string[j][i];
            } else {
                cur_fragment_is_missing = true;
                break;
            }
        }
        if (cur_fragment_is_missing) {
            sstream.str("");
            missing_profile_seq_index.push_back(i);
            continue;
        }
        file << sstream.str() << "\n";
        selected_profile_seq_index.push_back(i);
        sstream.str("");
    }
    file.close();
    // use the alignment tool to align the sequences.
    std::string cmd = "";
#if (defined(__linux__))
    cmd.append("ext/profile-two-align/linux/profile_two_align")
        .append(" -q ")
        .append(seq_file_name)
        .append(" -p ")
        .append(profile_file_name)
        .append(" -F")
        .append(" -D")
        .append(" 2> /dev/null");
#else
    cmd.append("ext\\profile-two-align\\win\\profile_two_align.exe")
        .append(" -q ")
        .append(seq_file_name)
        .append(" -p ")
        .append(profile_file_name)
        .append(" -F")
        .append(" -D")
        .append(" 2>NUL");
#endif

    int res = system(cmd.c_str());
    if (res != 0) {
#if DEBUG
        std::string out = "Warning: Seq-Profile alignment may result in errors and may produce invalid results.";
        print_table_line(out);
#endif
    }
    // Define vectors for storing aligned sequences and their names
    std::vector<std::string> align_res;
    std::vector<std::string> align_res_name;
    // Read data from a file and store it in the above vectors
    read_data(profile_file_name.c_str(), align_res, align_res_name, false);

    for (uint_t i = 0; i < selected_profile_seq_index.size(); i++) {
        concat_string[left_index][selected_profile_seq_index[i]] = align_res[i];
    }
    concat_string[left_index][seq_index] = align_res[align_res.size() - 1];
#if DEBUG
    std::cout << align_res[align_res.size() - 1] << std::endl;
#endif // DEBUG

    // Loop through missing profile sequence indices and update concat_string accordingly
    for (uint_t i = 0; i < missing_profile_seq_index.size(); i++) {
        for (uint_t j = left_index + 1; j <= right_index; j++) {
            concat_string[left_index][missing_profile_seq_index[i]] += concat_string[j][missing_profile_seq_index[i]];
        }
    }

    fragment_len[left_index] = concat_string[left_index][seq_index].length();
    // Update concat_range for all indices
    for (uint_t i = 0; i < data.size(); i++) {
        bool flag = false;
        for (uint_t j = left_index; j <= right_index; j++) {
            if (concat_range[j][i].first == -1) {
                flag = true;
            }
        }
        if (flag) {
            concat_range[left_index][i].first = -1;
            concat_range[left_index][i].second = -1;
            continue;
        }
        if (left_index > 0) {
            concat_range[left_index][i].first = concat_range[left_index - 1][i].first + concat_range[left_index - 1][i].second;
        } else {
            concat_range[left_index][i].first = 0;
        }
        if (right_index + 1 < concat_range.size()) {
            concat_range[left_index][i].second = concat_range[right_index + 1][i].first - concat_range[left_index][i].first;
        } else {
            concat_range[left_index][i].second = data[i].length() - concat_range[left_index][i].first;
        }
    }
    // Erase indices from concat_string, concat_range, and fragment_len vectors
    std::vector<std::vector<std::string>>::iterator str_it = concat_string.begin() + left_index + 1;
    std::vector<std::vector<std::pair<int_t, int_t>>>::iterator range_it = concat_range.begin() + left_index + 1;
    std::vector<uint_t>::iterator fragment_it = fragment_len.begin() + left_index + 1;
    int_t erase_len = right_index - left_index;
    while (erase_len > 0 && str_it != concat_string.end()) {
        str_it = concat_string.erase(str_it);
        range_it = concat_range.erase(range_it);
        fragment_it = fragment_len.erase(fragment_it);
        erase_len--;
    }
    // Remove files from disk
    if (remove(seq_file_name.c_str()) != 0) {
        std::cerr << "Error deleting file " << seq_file_name << std::endl;
    }
    if (remove(profile_file_name.c_str()) != 0) {
        std::cerr << "Error deleting file " << profile_file_name << std::endl;
    }
    return concat_string.begin() + left_index;
}

/**
 * Refines the given data by removing leading and trailing spaces.
 *
 * This function performs refinement on the given data, where the refinement process involves removing a certain
 * number of leading and trailing spaces (represented by '-') from the start and end of each string in the data. The
 * number of spaces to remove is determined by the string in the data with the smallest sum of leading and trailing
 * spaces.
 *
 * @param data1 A reference to the first vector of strings to refine.
 * @param data2 A reference to the second vector of strings to refine.
 */
void refinement(std::vector<std::string> &data1, std::vector<std::string> &data2) {
    // The number of spaces to remove from each string.
    int spaceToRemove = INT_MAX;

    // Determine the minimum number of spaces to remove.
    for (uint_t i = 0; i < data1.size(); ++i) {
        std::string &str1 = data1[i];
        std::string &str2 = data2[i];
        int spaceCount1 = 0, spaceCount2 = 0;

        // Count the number of trailing spaces in str1.
        while (str1.size() - 1 - spaceCount1 < str1.size() && str1[str1.size() - 1 - spaceCount1] == '-') {
            ++spaceCount1;
        }
        // Count the number of leading spaces in str2.
        while (static_cast<size_t>(spaceCount2) < str2.size() && str2[spaceCount2] == '-') {
            ++spaceCount2;
        }

        // The total number of leading and trailing spaces in str1 and str2.
        int totalSpace = spaceCount1 + spaceCount2;
        // Update the number of spaces to remove if necessary.
        spaceToRemove = spaceToRemove < totalSpace ? spaceToRemove : totalSpace;
    }

    // Perform the refinement process on each string in the data.
    for (uint_t i = 0; i < data1.size(); ++i) {
        std::string &str1 = data1[i];
        std::string &str2 = data2[i];
        int removedSpace = 0;

        // Remove trailing spaces from str1.
        while (str1.size() > 0 && removedSpace < spaceToRemove) {
            if (str1[str1.size() - 1] != '-')
                break;
            str1.pop_back();
            ++removedSpace;
        }
        // Remove leading spaces from str2.
        while (str2.size() > 0 && removedSpace < spaceToRemove) {
            if (str2[0] != '-')
                break;
            str2.erase(0, 1);
            ++removedSpace;
        }
    }
}

/**
 * @brief Concatenate two sets of sequence data (chain and parallel) into a single set of concatenated data.
 * @param chain_string A vector of vectors containing the chain sequence data.
 * @param parallel_string A vector of vectors containing the parallel sequence data.
 * @return std::vector<std::vectorstd::string> A vector of vectors containing the concatenated sequence data.
 */
std::vector<std::vector<std::string>> concat_chain_and_parallel(std::vector<std::vector<std::string>> &chain_string,
                                                                std::vector<std::vector<std::string>> &parallel_string) {
    // Determine the number of sequences in each set of data.
    uint_t seq_num = parallel_string[0].size();
    // Determine the number of sets of chain and parallel data.
    uint_t chain_num = chain_string.size();
    uint_t parallel_num = parallel_string.size();
    // Create a vector to hold the concatenated data.
    std::vector<std::vector<std::string>> concated_data(chain_num + parallel_num);
    uint_t count = 0;

    for (uint_t i = 0; i < parallel_num; ++i) {
        if (parallel_string[i].size() != seq_num) {
            // log + correção defensiva
            // (mesma lógica de resize mostrada acima)
        }
    }

    // Loop through each set of parallel data and add it to the concatenated data.
    for (uint_t i = 0; i < parallel_num; i++) {
        for (uint_t j = 0; j < seq_num; j++) {
            concated_data[count].push_back(parallel_string[i][j]);
            if (i < chain_num) {
                concated_data[count + 1].push_back(chain_string[i][j]);
            }
        }
        count += 2;
    }
    for (uint_t i = 0; i < chain_num + parallel_num - 1; i++) {
        refinement(concated_data[i], concated_data[i + 1]);
    }
    return concated_data;
}

/**
 * @brief Concatenates two sets of ranges in a chain and parallel manner
 * Given two vectors of vectors, chain and parallel, this function concatenates the ranges
 * in a chain and parallel manner. Specifically, it takes the i-th range from each vector in parallel
 * and concatenates them into a single vector. Then, it takes the i-th range from the chain vector
 * and concatenates it with the previous vector to form a new concatenated vector. The resulting
 * concatenated vector is stored in a new vector of vectors and returned.
 * @param chain A vector of vectors representing the chains of ranges to concatenate
 * @param parallel A vector of vectors representing the parallel ranges to concatenate
 * @return A new vector of vectors representing the concatenated ranges
 */
std::vector<std::vector<std::pair<int_t, int_t>>>
concat_chain_and_parallel_range(std::vector<std::vector<std::pair<int_t, int_t>>> &chain,
                                std::vector<std::vector<std::pair<int_t, int_t>>> &parallel) {
    uint_t seq_num = chain.size();
    uint_t chain_num = chain[0].size();
    uint_t parallel_num = parallel.size();
    std::vector<std::vector<std::pair<int_t, int_t>>> concated_range(chain_num + parallel_num);

    uint_t count = 0;
    for (uint_t i = 0; i < parallel_num; i++) {
        for (uint_t j = 0; j < seq_num; j++) {
            concated_range[count].push_back(parallel[i][j]);
            if (i < chain_num) {
                concated_range[count + 1].push_back(chain[j][i]);
            }
        }
        count += 2;
    }
    return concated_range;
}

/**
 * @brief Get the length of the first non-empty string in each row of a 2D vector of strings
 * @param concat_string The input 2D vector of strings
 * @return A vector of unsigned integers representing the length of the first non-empty string in each row
 */
std::vector<uint_t> get_first_nonzero_lengths(const std::vector<std::vector<std::string>> &concat_string) {
    std::vector<uint_t> result;
    for (const auto &row : concat_string) {
        uint_t length = 0;
        for (const auto &str : row) {
            if (!str.empty()) {
                length = str.length();
                break;
            }
        }
        result.push_back(length);
    }
    return result;
}
