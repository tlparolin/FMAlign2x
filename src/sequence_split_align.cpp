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
 * @brief Split sequences into chain-defined regions, perform parallel alignment,
 *        and produce a fully concatenated final multiple alignment.
 * This function orchestrates the entire FMAlign2 parallel alignment pipeline.
 * @param data  A vector of input sequences to be aligned.
 * @param name  Sequence identifiers in the same order as `data`.
 * @param chain A vector of chain (MEM-expanded) coordinate pairs, one vector
 *              per sequence, representing the initial anchors used for splitting.
 * @param pool  A ThreadPool instance used for all parallelizable steps.
 * @returns void (writes output to file specified in global_args).
 */
void split_and_parallel_align(std::vector<std::string> data, std::vector<std::string> name,
                              std::vector<std::vector<std::pair<int_t, int_t>>> chain, ThreadPool &pool) {
    if (global_args.verbose) {
        std::cout << "#                Parallel Aligning...                       #" << std::endl;
        print_table_divider();
    }

    Timer timer;
    uint_t chain_num = chain[0].size();

    // Holds the expanded (SW refined) chain regions for each block.
    std::vector<std::vector<std::string>> chain_string(chain_num);

    // -------------------------------------------------------------------------
    // 1. EXPAND CHAINS (local Smith–Waterman)
    // -------------------------------------------------------------------------
    std::vector<ExpandChainParams> params(chain_num);

    for (uint_t i = 0; i < chain_num; i++) {
        // Prepare parameter set for each expansion job
        params[i].data = &data;
        params[i].chain = &chain;
        params[i].chain_index = i;
        params[i].result_store = chain_string.begin() + i;
    }

    // Run expansion either in parallel or sequentially depending on coverage policy.
    if (global_args.min_seq_coverage == 1) {
        // Submit each expansion task to the pool
        for (uint_t i = 0; i < chain_num; i++) {
            pool.add_task([i, &params]() { expand_chain(&params[i]); });
        }
        pool.wait_for_tasks();
    } else {
        // Sequential fallback (used for partial-coverage strategies)
        for (uint_t i = 0; i < chain_num; i++) {
            expand_chain(&params[i]);
        }
    }

    params.clear();

    if (global_args.verbose) {
        print_table_line(std::format("SW expand time: {:.2f} seconds.", timer.elapsed_time()));
    }

    timer.reset();

    // -------------------------------------------------------------------------
    // 2. COMPUTE PARALLEL ALIGNMENT RANGES
    // -------------------------------------------------------------------------
    // These ranges define regions between chain anchors that must be aligned.
    std::vector<std::vector<std::pair<int_t, int_t>>> parallel_align_range = get_parallel_align_range(data, chain);

    // -------------------------------------------------------------------------
    // 3. PERFORM SPOA-BASED PARALLEL ALIGNMENTS
    // -------------------------------------------------------------------------
    // This step includes block splitting, SPOA calls, memory limits,
    // and multi-thread orchestration inside preprocess_parallel_blocks().
    std::vector<std::vector<std::string>> spoa_parallel_string = preprocess_parallel_blocks(data, parallel_align_range, pool);

    if (global_args.verbose) {
        print_table_line(std::format("Parallel Align Time: {:.2f} seconds.", timer.elapsed_time()));
    }

    timer.reset();

    // -------------------------------------------------------------------------
    // 4. CONCATENATE CHAIN + PARALLEL REGIONS
    // -------------------------------------------------------------------------
    // Builds full-length alignment in block form.
    auto concat_blocks = concat_chain_and_parallel(chain_string, spoa_parallel_string);

    // -------------------------------------------------------------------------
    // 5. FINAL OUTPUT (profile reconstruction & writing)
    // -------------------------------------------------------------------------
    concat_alignment_from_blocks(concat_blocks, name);

    if (global_args.verbose) {
        print_table_line(std::format("Seq-profile time: {:.2f} seconds.", timer.elapsed_time()));
        print_table_divider();
    }
}

static std::vector<std::string> merge_alignments(const std::vector<std::string> &msa1, const std::vector<std::string> &msa2) {
    // Caso trivial: MSA1 vazio
    if (msa1.empty() || msa1[0].empty()) {
        return msa2;
    }

    // Caso trivial: MSA2 vazio
    if (msa2.empty() || msa2[0].empty()) {
        return msa1;
    }

    const size_t len1 = msa1[0].size();
    const size_t len2 = msa2[0].size();

    // Se os comprimentos já coincidem → concatena direto
    if (len1 == len2) {
        std::vector<std::string> result;
        result.reserve(msa1.size() + msa2.size());
        result.insert(result.end(), msa1.begin(), msa1.end());
        result.insert(result.end(), msa2.begin(), msa2.end());
        return result;
    }

    // Determinar qual é maior e pad no outro
    const bool msa1_is_longer = (len1 > len2);
    const size_t target_len = msa1_is_longer ? len1 : len2;

    const std::vector<std::string> &base = msa1_is_longer ? msa1 : msa2;
    const std::vector<std::string> &pad_src = msa1_is_longer ? msa2 : msa1;

    std::vector<std::string> result;
    result.reserve(base.size() + pad_src.size());

    // Copia a MSA maior
    result.insert(result.end(), base.begin(), base.end());

    // Adiciona a menor com pad
    for (const auto &seq : pad_src) {
        std::string padded = seq;
        padded.resize(target_len, '-');
        result.push_back(std::move(padded));
    }

    return result;
}

bool will_use_clustering(size_t num_sequences, size_t cluster_size) {
    return num_sequences > cluster_size * 2; // Se > 2x o tamanho do cluster
}

std::vector<std::string> align_smart(const std::vector<std::string> &sequences, size_t cluster_size) {

    // Caso 1: Vazio
    if (sequences.empty()) {
        return {};
    }

    // Caso 2: Uma sequência
    if (sequences.size() == 1) {
        return sequences;
    }

    // Caso 3: Poucas sequências - usar SPOA direto
    if (!will_use_clustering(sequences.size(), cluster_size)) {
        return run_spoa_local(sequences);
    }

    // Caso 4: MUITAS sequências - usar clustering
    // Dividir sequências em clusters
    std::vector<std::vector<std::string>> clusters;
    for (size_t i = 0; i < sequences.size(); i += cluster_size) {
        size_t end = std::min(i + cluster_size, sequences.size());
        std::vector<std::string> cluster(sequences.begin() + i, sequences.begin() + end);
        clusters.push_back(cluster);
    }

    // Alinhar cada cluster com SPOA
    std::vector<std::vector<std::string>> aligned_clusters;
    for (size_t i = 0; i < clusters.size(); ++i) {
        auto aligned = run_spoa_local(clusters[i]);
        aligned_clusters.push_back(aligned);
    }

    // Merge: usar o primeiro cluster como base, adicionar os outros
    std::vector<std::string> result = aligned_clusters[0];

    for (size_t i = 1; i < aligned_clusters.size(); ++i) {
        result = merge_alignments(result, aligned_clusters[i]);
    }
    return result;
}

/**
 * @brief Concatenate aligned blocks into a final FASTA-formatted alignment output.
 *
 * This function takes a vector of aligned blocks — where each block represents
 * a multiple sequence alignment (MSA) of a MEM region or an inter-MEM interval —
 * and reconstructs the full-length alignment by concatenating these blocks in order.
 *
 * After concatenation, the function validates that all resulting sequences have
 * identical lengths (as expected for a valid MSA) and writes the final alignment
 * in FASTA format to the output file defined in global settings.
 *
 * @param blocks A vector of aligned blocks. Each block contains one aligned string
 *        per input sequence.
 * @param names The names of the sequences in the final alignment.
 * @return void (writes the aligned sequences to disk)
 *
 * @throws std::runtime_error If validation fails or the output file cannot be opened.
 */
void concat_alignment_from_blocks(const std::vector<std::vector<std::string>> &blocks, const std::vector<std::string> &names) {
    // =========================================================================
    // Input Validation
    // =========================================================================
    if (blocks.empty()) {
        // Nothing to concatenate; likely an upstream failure or empty input
        if (global_args.verbose) {
            print_table_line("Warning: Empty blocks provided for concatenation");
        }
        return;
    }

    size_t seq_num = names.size();

    if (seq_num == 0) {
        throw std::runtime_error("Error: No sequence names provided.");
    }

    // Ensure each block contains one entry per sequence
    for (size_t b = 0; b < blocks.size(); ++b) {
        if (blocks[b].size() != seq_num) {
            throw std::runtime_error(std::format("Error: Block {} has {} sequences, but {} were expected.", b, blocks[b].size(), seq_num));
        }
    }

    // =========================================================================
    // Concatenate Blocks
    // =========================================================================
    // Pre-allocate output strings to avoid repeated reallocations
    std::vector<std::string> final_alignment(seq_num);

    // Estimate total alignment length for memory reservation
    size_t total_length = 0;
    for (const auto &block : blocks) {
        for (const auto &seq : block) {
            total_length += seq.size();
        }
    }

    // Reserve space for each final sequence
    for (size_t s = 0; s < seq_num; ++s) {
        final_alignment[s].reserve(total_length / seq_num + 100);
    }

    // Concatenate aligned blocks in order
    for (const auto &block : blocks) {
        for (size_t s = 0; s < seq_num; ++s) {
            // Avoid adding empty segments
            if (!block[s].empty()) {
                final_alignment[s] += block[s];
            }
        }
    }

    // =========================================================================
    // Validation After Concatenation
    // =========================================================================
    // All final sequences must have equal length (property of an MSA)
    size_t expected_len = final_alignment[0].size();

    for (size_t s = 1; s < seq_num; ++s) {
        if (final_alignment[s].size() != expected_len) {
            throw std::runtime_error(std::format("Error: After concatenation, sequence {} has length {} while "
                                                 "sequence 0 has length {}.",
                                                 s, final_alignment[s].size(), expected_len));
        }
    }

    // =========================================================================
    // Write FASTA Output
    // =========================================================================
    std::string output_path = global_args.output_path;
    std::ofstream output_file(output_path);

    if (!output_file.is_open()) {
        throw std::runtime_error(std::format("Error: Cannot open output file '{}'.", output_path));
    }

    // Write each sequence in FASTA format
    for (size_t s = 0; s < seq_num; ++s) {
        output_file << ">" << names[s] << "\n" << final_alignment[s] << "\n";
    }
    // output_file will close automatically when going out of scope
}

/**
 * @brief Extract sequence fragments based on per-sequence ranges.
 *
 * This function receives the original sequences and a vector of (start, length)
 * ranges, and extracts one fragment from each sequence using these coordinates.
 * Invalid ranges (out of bounds, negative values, or empty segments) produce
 * empty fragments. The function also returns the maximum fragment length extracted.
 *
 * @param data The original input sequences.
 * @param range A vector of (start, length) pairs defining the extraction window
 *        for each sequence.
 * @param fragments Output vector where extracted fragments will be stored.
 * @return The maximum length among all extracted fragments.
 */
size_t extract_fragments(const std::vector<std::string> &data, const std::vector<std::pair<int, int>> &range,
                         std::vector<std::string> &fragments) {
    size_t max_len = 0;
    fragments.resize(data.size());

    for (uint_t s = 0; s < data.size(); ++s) {

        // Validate that a range exists for this sequence index
        if (s >= range.size()) {
            fragments[s].clear();
            continue;
        }

        auto [start, len] = range[s];

        // Validate start and length values
        if (start < 0 || len <= 0 || start >= static_cast<int>(data[s].size())) {
            fragments[s].clear();
            continue;
        }

        // Clamp length to avoid going past sequence boundaries
        if (start + len > static_cast<int>(data[s].size())) {
            len = data[s].size() - start;
        }

        // Extract the fragment
        fragments[s] = data[s].substr(start, len);

        // Track maximum fragment length
        max_len = std::max(max_len, fragments[s].size());
    }

    return max_len;
}

/**
 * @brief Check whether all extracted fragments are identical.
 *
 * This function verifies if all fragments in the input vector are equal
 * to the first non-empty fragment. If the fragment list is empty or the
 * first fragment is empty, the function returns false, since equality
 * cannot be meaningfully evaluated.
 *
 * @param fragments A vector of sequence fragments.
 * @return true if all fragments are identical and non-empty, false otherwise.
 */
bool all_fragments_equal(const std::vector<std::string> &fragments) {
    // If the list is empty or the first fragment is empty,
    // we cannot meaningfully compare equality.
    if (fragments.empty() || fragments[0].empty())
        return false;

    const std::string &first = fragments[0];

    // Check if every fragment matches the first one exactly
    return std::all_of(fragments.begin(), fragments.end(), [&first](const std::string &s) { return s == first; });
}

/**
 * @brief Extract a sub-fragment from an existing fragment.
 *
 * This function extracts a substring of `fragment` between the positions
 * [pos, end). If the indices are invalid (out of bounds, inverted range,
 * or empty input), an empty string is returned. Bounds are clamped to the
 * fragment size to ensure safe extraction.
 *
 * @param fragment The original fragment from which to extract.
 * @param pos The starting position (inclusive).
 * @param end The ending position (exclusive).
 * @return A substring of the fragment within the specified range, or an empty
 *         string if the range is invalid.
 */
std::string extract_sub_fragment(const std::string &fragment, size_t pos, size_t end) {
    // Validate input fragment and ensure start position is within bounds
    if (fragment.empty() || pos >= fragment.size()) {
        return "";
    }

    // Ensure valid interval (pos must be < end)
    if (pos >= end) {
        return "";
    }

    // Clamp the start and end to the fragment boundaries
    size_t real_start = std::min(pos, fragment.size());
    size_t real_end = std::min(end, fragment.size());

    // Extract substring only if the resulting range is valid
    if (real_start < real_end) {
        return fragment.substr(real_start, real_end - real_start);
    }

    return "";
}

/**
 * @brief Merges aligned sub-blocks into a final full-block alignment by trimming overlaps.
 *
 * This function performs a **non-alignment merge** of multiple subdivision MSAs.
 * It assumes that each sub-block represents a correctly aligned portion of the
 * original sequences with an intentional overlap region. The merge strategy:
 *
 *   1. Copy the first aligned block entirely.
 *   2. For each subsequent block, remove `overlap` columns from the end of the
 *      already merged alignment.
 *   3. Append the next aligned block as-is.
 *
 * No re-alignment is performed. Only string trimming and concatenation.
 *
 * @param sub_results A vector of MSAs (each MSA is a vector<string> of equal size).
 * @param original The original unaligned fragments (used only as fallback reference).
 * @param overlap The number of alignment columns to trim when merging consecutive blocks.
 * @return A vector<string> containing the merged full-block alignment.
 *
 * @note This function assumes that each sub-block has already been aligned
 *       independently (e.g., using SPOA). The merging is purely structural.
 */
std::vector<std::string> merge_subdivisions_simple(const std::vector<std::vector<std::string>> &sub_results,
                                                   const std::vector<std::string> &original, size_t overlap) {

    size_t seq_num = original.size();
    std::vector<std::string> merged(seq_num, "");

    // --- Validation: empty input → return original fragments ---
    if (sub_results.empty()) {
        return original;
    }

    // --- Validation: first block must match expected number of sequences ---
    if (sub_results[0].size() != seq_num) {
        // Corrupted input; safest recovery is to return original
        return original;
    }

    // Copy the first aligned sub-block entirely
    merged = sub_results[0];

    // Merge remaining sub-blocks
    for (size_t sub_idx = 1; sub_idx < sub_results.size(); ++sub_idx) {

        // Validate block size before using it
        if (sub_results[sub_idx].size() != seq_num) {
            // Skip corrupted block; do not fail the entire merge
            continue;
        }

        const auto &sub = sub_results[sub_idx];

        // Merge each sequence independently
        for (size_t seq = 0; seq < seq_num; ++seq) {

            // Bounds check (stronger than needed but safe)
            if (seq >= merged.size() || seq >= sub.size()) {
                continue;
            }

            const std::string &prev = merged[seq];
            const std::string &curr = sub[seq];

            // If current block is empty, nothing to merge
            if (curr.empty()) {
                continue;
            }

            // Compute how much overlap to trim
            size_t overlap_len = std::min(overlap, prev.size());

            // Trim overlap region from the end of the already-merged alignment
            if (overlap_len > 0 && prev.size() >= overlap_len) {
                merged[seq].erase(prev.size() - overlap_len);
            }

            // Append the current sub-block
            merged[seq] += curr;
        }
    }

    return merged;
}

/**
 * @brief Preprocess alignment blocks using SPOA, with optional subdivision for large regions.
 *
 * This function performs a fast pre-alignment step for each block defined in
 * `parallel_align_range`. The processing strategy is:
 *
 *   1. **Identical fragments** → direct copy (no alignment needed)
 *   2. **Small blocks (< global_args.max_block_size)** → single-step SPOA alignment
 *   3. **Large blocks (>= global_args.max_block_size)** → subdivide, align sub-blocks independently
 *      and merge using overlap trimming
 *
 * This pre-alignment reduces load on heavier external aligners (e.g., MAFFT, HAlign)
 * and is optimized for ultralong genomic sequences between MEM anchors.
 *
 * All blocks are scheduled immediately to the provided `ThreadPool` so that thread
 * utilization remains high even when individual alignment tasks vary in cost.
 *
 * @param data Original input sequences.
 * @param parallel_align_range A list of (start, length) ranges for each block and each sequence.
 * @param pool Thread pool used to run SPOA tasks concurrently.
 * @return A list of aligned blocks (no external reprocessing required).
 *
 * @note Subdivision overlap is fixed to global_args.overlap_size to maintain block continuity.
 * @note SPOA failures are caught, and original fragments are used as fallback.
 */
std::vector<std::vector<std::string>> preprocess_parallel_blocks(const std::vector<std::string> &data,
                                                                 const std::vector<std::vector<std::pair<int, int>>> &parallel_align_range,
                                                                 ThreadPool &pool) {
    uint_t parallel_num = parallel_align_range.size();
    uint_t seq_num = data.size();

    // Output vector: one aligned block per task
    std::vector<std::vector<std::string>> result(parallel_num, std::vector<std::string>(seq_num));

    // Statistics for diagnostics
    std::atomic<uint_t> count_exact(0);
    std::atomic<uint_t> count_spoa_direct(0);
    std::atomic<uint_t> count_spoa_subdivided(0);

    // Capture pointers to avoid expensive capturing in lambda
    const auto *data_ptr = &data;
    const auto *range_ptr = &parallel_align_range;

    // ========================================================================
    // PASS 1 — Launch ALL block tasks immediately
    // ========================================================================

    for (uint_t block_idx = 0; block_idx < parallel_num; ++block_idx) {

        pool.add_task([block_idx, data_ptr, range_ptr, &result, &count_exact, &count_spoa_direct, &count_spoa_subdivided]() {
            const auto &data_ref = *data_ptr;
            const auto &range_ref = *range_ptr;

            // Validate index (safety check)
            if (block_idx >= range_ref.size()) {
                return;
            }

            const auto &range = range_ref[block_idx];

            // Extract block fragments (per-sequence)
            std::vector<std::string> fragments(data_ref.size());
            size_t max_len = extract_fragments(data_ref, range, fragments);

            // ================================================================
            // CASE 1: All fragments identical → direct copy
            // ================================================================
            if (all_fragments_equal(fragments)) {
                std::fill(result[block_idx].begin(), result[block_idx].end(), fragments[0]);
                count_exact.fetch_add(1, std::memory_order_relaxed);
                return;
            }

            // ================================================================
            // CASE 2: Small block (< global_args.max_block_size) → direct SPOA alignment
            // ================================================================
            if (max_len <= global_args.max_block_size) {
                result[block_idx] = align_smart(fragments, 500);
                count_spoa_direct.fetch_add(1, std::memory_order_relaxed);
                return;
            }

            // ================================================================
            // CASE 3: Large block → subdivide + run SPOA on each sub-block
            // ================================================================
            count_spoa_subdivided.fetch_add(1, std::memory_order_relaxed);

            size_t stride = std::max(global_args.max_block_size - global_args.overlap_size, global_args.max_block_size / 2);
            size_t num_subblocks = (max_len + stride - 1) / stride;

            // If something went wrong and no subdivisions computed → fallback
            if (num_subblocks == 0) {
                result[block_idx] = fragments;
                return;
            }

            std::vector<std::vector<std::string>> sub_results(num_subblocks);

            size_t sub_idx = 0;

            // Process each subdivision in order
            for (size_t pos = 0; pos < max_len; pos += stride) {

                size_t end = std::min(pos + global_args.max_block_size + global_args.overlap_size, max_len);

                // Extract sub-fragments for this subdivision
                std::vector<std::string> sub_frags(data_ref.size());
                for (uint_t s = 0; s < data_ref.size(); ++s) {
                    sub_frags[s] = extract_sub_fragment(fragments[s], pos, end);
                }

                // Try to align subdivision with SPOA
                try {
                    // sub_results[sub_idx] = run_spoa_local(sub_frags);
                    size_t cluster_size = sub_frags.size() > 2000 ? 500 : 200;
                    sub_results[sub_idx] = align_smart(sub_frags, cluster_size);
                } catch (...) {
                    // SPOA failure → store raw (unaligned) fragments
                    sub_results[sub_idx] = sub_frags;
                }

                if (end >= max_len) {
                    break;
                }

                ++sub_idx;
            }

            // Merge all aligned subdivisions into a single block
            result[block_idx] = merge_subdivisions_simple(sub_results, fragments, global_args.overlap_size);
        });
    }

    // ========================================================================
    // PASS 2 — Wait for ALL job completions
    // ========================================================================
    pool.wait_for_tasks();

    // Output diagnostic table (optional)
    if (global_args.verbose) {
        print_table_line("Exact-copy blocks:      " + std::to_string(count_exact.load()));
        print_table_line("Direct SPOA blocks:     " + std::to_string(count_spoa_direct.load()));
        print_table_line("Subdivided SPOA blocks: " + std::to_string(count_spoa_subdivided.load()));
    }

    return result;
}

/**
 * @brief Run a SPOA (Simd Partial Order Alignment) multiple-sequence alignment locally.
 *
 * This helper function performs MSA on a set of sequences, typically used for
 * subdivided alignment blocks. Empty sequences are allowed and handled: they
 * are ignored during graph construction and later reinserted as gap-only
 * sequences in the final MSA.
 *
 * If SPOA fails at any point, the function safely returns the input sequences
 * unchanged, ensuring robustness in the alignment pipeline.
 *
 * @param seqs Input sequences to align (may include empty strings).
 * @return A vector of aligned sequences of equal length; empty input yields an empty result.
 */
std::vector<std::string> run_spoa_local(const std::vector<std::string> &seqs) {

    // If no sequences were provided, return empty
    if (seqs.empty())
        return {};

    try {
        // Create SPOA alignment engine (global parameters)
        auto aligner = spoa::AlignmentEngine::Create(spoa::AlignmentType::kNW, 5, -4, -8);

        spoa::Graph graph;
        std::vector<size_t> non_empty_indices;
        non_empty_indices.reserve(seqs.size());

        // Add only non-empty sequences into the SPOA graph
        for (size_t i = 0; i < seqs.size(); ++i) {
            if (seqs[i].empty())
                continue;

            try {
                auto alignment = aligner->Align(seqs[i], graph);
                graph.AddAlignment(alignment, seqs[i]);
                non_empty_indices.push_back(i);
            } catch (...) {
                // If a sequence causes an exception during alignment, skip it
                continue;
            }
        }

        // If all sequences were empty, return a vector of empty fragments
        if (non_empty_indices.empty()) {
            return std::vector<std::string>(seqs.size(), "");
        }

        // Generate MSA from SPOA graph
        auto msa_compact = graph.GenerateMultipleSequenceAlignment();
        size_t aln_len = msa_compact.empty() ? 0 : msa_compact[0].size();

        // Validate generated MSA
        if (aln_len == 0) {
            // If SPOA produced no alignment, return original
            return seqs;
        }

        // Expand the MSA back to the original number of sequences
        // Empty original sequences become full-gap strings
        std::vector<std::string> msa(seqs.size(), std::string(aln_len, '-'));

        for (size_t k = 0; k < non_empty_indices.size(); ++k) {
            if (k < msa_compact.size()) {
                msa[non_empty_indices[k]] = std::move(msa_compact[k]);
            }
        }

        return msa;

    } catch (...) {
        // If SPOA fails entirely, return original sequences
        return seqs;
    }
}

/**
 * @brief Expands a MEM chain block for all sequences using local Smith–Waterman alignment.
 *
 * This function retrieves the chain block at the specified chain index and attempts
 * to reconstruct the missing aligned fragment for each sequence. The first available
 * non-empty chain entry is used as the query sequence (all other sequences attempt
 * to align against it).
 *
 * Alignment strategy:
 * - If the current sequence has a valid chain interval → copy the corresponding fragment.
 * - If the interval is missing → determine a reference interval based on neighboring
 *   chains, extract that reference fragment, and run Smith–Waterman to align the query
 *   against it.
 *
 * Neighbor-based reference window:
 * - Left boundary  = end of previous valid chain block (or 0 if none)
 * - Right boundary = start of next valid chain block (or end of sequence if none)
 *
 * The alignment result updates:
 * - `aligned_fragment[i]` with the aligned sequence fragment
 * - `chain[i][chain_index]` with the new alignment coordinates when alignment succeeds
 *
 * @param arg Pointer to an ExpandChainParams struct (opaque void* to fit pthread signature).
 * @return Always returns nullptr.
 */
void *expand_chain(void *arg) {
    // Cast raw pointer to parameter struct
    auto *ptr = static_cast<ExpandChainParams *>(arg);

    // References to external structures (no copies)
    const auto &data = *(ptr->data);
    auto &chain = *(ptr->chain);

    const size_t chain_index = ptr->chain_index;
    const size_t seq_num = data.size();
    const size_t chain_num = chain[0].size();

    // Smith–Waterman alignment engine (Striped SIMD implementation)
    StripedSmithWaterman::Aligner aligner;
    StripedSmithWaterman::Filter filter;
    StripedSmithWaterman::Alignment alignment;

    // -------------------------------------------------------------------------
    // STEP 1 — Determine the QUERY sequence
    // -------------------------------------------------------------------------

    std::string query;
    int_t query_length = 0;
    bool found_query = false;

    // The first sequence with a valid MEM interval becomes the query
    for (size_t i = 0; i < seq_num; ++i) {
        if (chain[i][chain_index].first != -1) {
            query_length = chain[i][chain_index].second;
            query = data[i].substr(chain[i][chain_index].first, query_length);
            found_query = true;
            break;
        }
    }

    // If no sequence contains a valid MEM: return empty alignment
    if (!found_query) {
        *(ptr->result_store) = std::vector<std::string>(seq_num, "");
        return nullptr;
    }

    // -------------------------------------------------------------------------
    // STEP 2 — Run alignment for all sequences
    // -------------------------------------------------------------------------

    std::vector<std::string> aligned_fragment(seq_num);

    for (size_t i = 0; i < seq_num; ++i) {

        int_t begin_pos = chain[i][chain_index].first;

        // =====================================================================
        // CASE A: Missing MEM for this sequence → must align using SW
        // =====================================================================
        if (begin_pos == -1) {

            size_t tmp_index = chain_index;
            // Mask length helps focus SW on local gaps (Smith–Waterman heuristic)
            int_t maskLen = std::max(query_length / 2, 15);

            // -----------------------------------------
            // Determine reference interval boundaries
            // -----------------------------------------
            size_t ref_begin_pos = 0;

            // If previous chain is valid → use its end as the reference start
            if (tmp_index > 0 && chain[i][tmp_index - 1].first != -1) {
                ref_begin_pos = chain[i][tmp_index - 1].first + chain[i][tmp_index - 1].second;
            }

            // If next chain is valid → use its start as the reference end
            size_t ref_end_pos = data[i].length() - 1;
            if (tmp_index < chain_num - 1 && chain[i][tmp_index + 1].first != -1) {
                ref_end_pos = chain[i][tmp_index + 1].first;
            }

            // Extract reference subsequence
            std::string ref = data[i].substr(ref_begin_pos, ref_end_pos - ref_begin_pos);

            // -----------------------------------------
            // Perform Smith–Waterman alignment
            // -----------------------------------------
            aligner.Align(query.c_str(), ref.c_str(), ref.size(), filter, &alignment, maskLen);

            // Store and interpret results
            auto p = store_sw_alignment(alignment, ref, query, aligned_fragment, i);

            // Update chain coordinates if alignment returned a valid match
            if (p.first != -1) {
                p.first += ref_begin_pos;  // Convert relative offset to global
                chain[i][chain_index] = p; // Save in chain structure
            }

        } else {
            // =================================================================
            // CASE B: Already has a valid MEM → simply copy the fragment
            // =================================================================
            aligned_fragment[i] = query;
        }
    }

    // -------------------------------------------------------------------------
    // STEP 3 — Store results and finish
    // -------------------------------------------------------------------------
    *(ptr->result_store) = std::move(aligned_fragment);
    return nullptr;
}

/**
 * @brief Convert a StripedSmithWaterman alignment into an aligned reference string and
 *        propagate required gaps into all previously aligned sequences.
 *
 * This function interprets the CIGAR operations generated by the SSW alignment and
 * reconstructs the alignment in a left-anchored fashion:
 *
 * - The aligned reference sequence is constructed into `res_store[seq_index]`.
 * - If a deletion in the query (`D`) occurs, gaps are inserted into the query and also
 *   retroactively into all *previous* aligned sequences to maintain column consistency.
 * - Soft clips (`S`) at the beginning or end extend the alignment window and may shift
 *   reference coordinates.
 *
 * Alignment failure conditions:
 * - If alignment.ref_begin < 0 → return (-1, -1)
 * - If soft-clipped region exceeds 80% of the alignment → considered unstable → fail
 *
 * The function returns:
 * - (ref_start, ref_length) where ref_start and ref_length correspond to the adjusted
 *   region of the reference covered by the alignment.
 *
 * @param alignment  Alignment result produced by StripedSmithWaterman.
 * @param ref        The reference sequence.
 * @param query      The query sequence (will be modified if gaps must be inserted).
 * @param res_store  Vector of aligned sequences; the aligned result is stored at seq_index.
 * @param seq_index  Index of the current sequence in the alignment block.
 * @return std::pair<int_t,int_t>  {ref_start, ref_length} or {-1,-1} on failure.
 */
std::pair<int_t, int_t> store_sw_alignment(StripedSmithWaterman::Alignment alignment, std::string &ref, std::string &query,
                                           std::vector<std::string> &res_store, uint_t seq_index) {
    // Reference to the CIGAR operations (SSW encodes them as 32-bit ints)
    const std::vector<unsigned int> &cigar = alignment.cigar;

    int_t ref_begin = alignment.ref_begin;
    int_t ref_end = alignment.ref_end;
    uint_t query_begin = 0; // query always begins at offset 0 in this implementation

    // -------------------------------------------------------------------------
    // FAILURE CASE: alignment not found
    // -------------------------------------------------------------------------
    if (ref_begin < 0) {
        res_store[seq_index].clear();
        return {-1, -1};
    }

    // -------------------------------------------------------------------------
    // Handle soft-clips and compute extended boundaries
    // -------------------------------------------------------------------------
    int_t S_count = 0;
    int_t total_length = alignment.query_end;
    int_t new_ref_begin = ref_begin;
    uint_t new_ref_end = ref_end;

    auto cigar_len = [](unsigned int c) { return cigar_int_to_len(c); };
    auto cigar_op = [](unsigned int c) { return cigar_int_to_op(c); };

    // Soft-clip at the beginning
    if (cigar_op(cigar.front()) == 'S') {
        int_t len = cigar_len(cigar.front());
        S_count += len;
        new_ref_begin = std::max<int_t>(0, new_ref_begin - len);
    }

    // Soft-clip at the end
    if (cigar_op(cigar.back()) == 'S') {
        int_t len = cigar_len(cigar.back());
        S_count += len;
        total_length += len;
        new_ref_end = std::min<uint_t>(ref.length() - 1, new_ref_end + len);
    }

    // Excessive soft-clipping → alignment considered unreliable
    if (S_count > static_cast<int_t>(std::ceil(0.8 * total_length))) {
        res_store[seq_index].clear();
        return {-1, -1};
    }

    // -------------------------------------------------------------------------
    // Reconstruct aligned reference sequence
    // -------------------------------------------------------------------------
    std::string aligned_result;
    uint_t p_ref = 0;   // position into reference alignment window
    uint_t p_query = 0; // position into query (including gaps)

    // -------------------------------------------------------------------------
    // Iterate through CIGAR operations and reconstruct alignment
    // -------------------------------------------------------------------------
    for (const auto &c : cigar) {

        char op = cigar_op(c);
        int_t len = cigar_len(c);

        switch (op) {

        // ================================================================
        // SOFT CLIP (S) — sequence outside alignment but may influence
        // reference anchoring and aligned output.
        // ================================================================
        case 'S': {
            if (&c == &cigar.front()) {
                // At start: extend backward
                int_t tmp = len;

                // If soft-clip reaches past the alignment start, pad with '-'
                if (ref_begin <= len) {
                    while (tmp > ref_begin) {
                        aligned_result += '-';
                        --tmp;
                    }
                }

                // Append clipped reference bases before alignment start
                for (int_t j = ref_begin - tmp; j < ref_begin; ++j) {
                    aligned_result += ref[j];
                }
            } else {
                // At end: append clipped region after ref_end
                int_t tmp = len;

                for (uint_t j = ref_end + 1; j < ref.length() && tmp > 0; ++j) {
                    aligned_result += ref[j];
                    --tmp;
                }

                // If ref has no more characters, pad remaining S with gaps
                while (tmp-- > 0)
                    aligned_result += '-';
            }

            p_query += len;
            break;
        }

        // ================================================================
        // MATCH / MISMATCH / EQUAL
        // ================================================================
        case 'M':
        case 'X':
        case '=': {
            aligned_result.append(ref.begin() + ref_begin + p_ref, ref.begin() + ref_begin + p_ref + len);
            p_ref += len;
            p_query += len;
            break;
        }

        // ================================================================
        // INSERTION (I) — insertion relative to reference → gap in reference
        // ================================================================
        case 'I': {
            aligned_result.append(len, '-'); // query has bases; ref has '-'
            p_query += len;
            break;
        }

        // ================================================================
        // DELETION (D) — deletion in query → reference bases appear here,
        // but a gap must be inserted in the query and *propagated to all
        // previous aligned sequences* to maintain global column consistency.
        // ================================================================
        case 'D': {
            // Append the deleted reference region directly
            aligned_result.append(ref.begin() + ref_begin + p_ref, ref.begin() + ref_begin + p_ref + len);

            p_ref += len;

            // Insert identical gaps into the query and all sequences before it
            std::string gaps(len, '-');

            // Insert in query
            query.insert(query_begin + p_query, gaps);

            // Propagate deletion to previously aligned sequences
            for (uint_t j = 0; j < seq_index; ++j) {
                if (!res_store[j].empty()) {
                    res_store[j].insert(query_begin + p_query, gaps);
                }
            }

            p_query += len;
            break;
        }
        }
    }

    // -------------------------------------------------------------------------
    // Store final aligned reference result
    // -------------------------------------------------------------------------
    res_store[seq_index] = std::move(aligned_result);

    // -------------------------------------------------------------------------
    // Return updated reference mapping
    // -------------------------------------------------------------------------
    return {new_ref_begin, new_ref_end - new_ref_begin + 1};
}

/**
 * @brief Compute the alignment ranges for each sequence between chain blocks.
 *
 * This function takes the detected chain blocks (MEM-expanded regions) for each
 * sequence and computes the ranges of "unaligned" segments that lie between these
 * chain blocks. The resulting structure is transposed so that each entry corresponds
 * to a block of parallel alignment, containing the ranges for all sequences.
 *
 * @param data  Vector of sequences to align.
 * @param chain Per-sequence chain blocks represented as pairs {begin, length}.
 * @return A transposed vector where each element is a block and contains the ranges
 *         {start, length} for all sequences in that block.
 */
std::vector<std::vector<std::pair<int_t, int_t>>> get_parallel_align_range(const std::vector<std::string> &data,
                                                                           const std::vector<std::vector<std::pair<int_t, int_t>>> &chain) {
    const uint_t seq_num = data.size();

    // Determine the maximum number of chain blocks across all sequences
    const uint_t chain_num =
        std::max_element(chain.begin(), chain.end(), [](const auto &a, const auto &b) { return a.size() < b.size(); })->size();

    // Preallocate memory for per-sequence ranges
    std::vector<std::vector<std::pair<int_t, int_t>>> parallel_align_range(seq_num);

    // Compute the per-sequence ranges
    for (uint_t i = 0; i < seq_num; i++) {
        const std::string &seq = data[i];
        const auto &seq_chain = chain[i];

        int_t last_pos = 0;
        auto &tmp_range = parallel_align_range[i];
        tmp_range.reserve(chain_num + 1);

        for (uint_t j = 0; j < chain_num; j++) {
            // Get chain entry or set to {-1,-1} if it does not exist
            std::pair<int_t, int_t> entry = (j < seq_chain.size()) ? seq_chain[j] : std::pair<int_t, int_t>{-1, -1};

            const int_t begin_pos = entry.first;

            if (begin_pos == -1) {
                tmp_range.emplace_back(-1, -1);
                last_pos = -1; // entering missing-state
            } else {
                if (last_pos == -1) {
                    // If previous block was missing, this parallel segment must also be missing
                    tmp_range.emplace_back(-1, -1);
                } else {
                    // Normal gap between blocks
                    tmp_range.emplace_back(last_pos, begin_pos - last_pos);
                }

                // Update last_pos safely
                last_pos = (entry.second >= 0) ? begin_pos + entry.second : begin_pos;
            }
        }

        // Handle trailing region after last chain block
        if (last_pos == -1) {
            tmp_range.emplace_back(-1, -1);
        } else {
            int_t remaining = static_cast<int_t>(seq.length()) - last_pos;
            tmp_range.emplace_back(last_pos, std::max<int_t>(remaining, 0));
        }
    }

    // Transpose result: block-major format
    std::vector<std::vector<std::pair<int_t, int_t>>> transpose_res;
    transpose_res.reserve(chain_num + 1);

    for (uint_t j = 0; j <= chain_num; j++) {
        std::vector<std::pair<int_t, int_t>> row;
        row.reserve(seq_num);

        for (uint_t i = 0; i < seq_num; i++) {
            row.push_back(j < parallel_align_range[i].size() ? parallel_align_range[i][j] : std::pair<int_t, int_t>{-1, -1});
        }

        transpose_res.push_back(std::move(row));
    }

    return transpose_res;
}

/**
 * @brief Concatenate chain-aligned blocks and parallel-aligned blocks into a single structure.
 *
 * This function interleaves parallel alignment blocks and chain blocks. For each index `i`,
 * the order becomes:
 *
 *     [ parallel_string[i], chain_string[i], parallel_string[i+1], chain_string[i+1], ... ]
 *
 * If one of the structures (parallel or chain) contains more blocks than the other,
 * the remaining blocks are appended at the end.
 *
 * @param chain_string    Vector of per-block sequence sets aligned by chain/MEM expansion.
 * @param parallel_string Vector of per-block sequence sets aligned in parallel.
 * @return A vector of vector<string> where each element represents a final concatenated block.
 */
std::vector<std::vector<std::string>> concat_chain_and_parallel(const std::vector<std::vector<std::string>> &chain_string,
                                                                const std::vector<std::vector<std::string>> &parallel_string) {
    // If both inputs are empty, return empty
    if (chain_string.empty() && parallel_string.empty())
        return {};

    // Determine number of sequences per block
    uint_t seq_num = 0;
    if (!parallel_string.empty())
        seq_num = parallel_string[0].size();
    else
        seq_num = chain_string[0].size();

    const uint_t chain_num = chain_string.size();
    const uint_t parallel_num = parallel_string.size();

    // Total number of resulting blocks = sum of both (interleaving)
    const uint_t total_blocks = chain_num + parallel_num;

    std::vector<std::vector<std::string>> concatenated;
    concatenated.resize(total_blocks);

    uint_t out_idx = 0; // index in output array

    const uint_t n = std::max(chain_num, parallel_num);

    // Interleave blocks: parallel[i], chain[i], ...
    for (uint_t i = 0; i < n; i++) {
        // Add parallel block if exists
        if (i < parallel_num) {
            auto &out = concatenated[out_idx++];
            out.reserve(seq_num);
            for (uint_t s = 0; s < seq_num; s++)
                out.push_back(parallel_string[i][s]);
        }

        // Add chain block if exists
        if (i < chain_num) {
            auto &out = concatenated[out_idx++];
            out.reserve(seq_num);
            for (uint_t s = 0; s < seq_num; s++)
                out.push_back(chain_string[i][s]);
        }
    }

    return concatenated;
}
