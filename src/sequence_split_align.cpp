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

#include "sequence_split_align.h"

/**
 * @brief Generates a random string of the specified length.
 * This function generates a random string of the specified length. The generated string
 * consists of lowercase English letters ('a' to 'z') for Linux platforms, and random bytes
 * for Windows platforms.
 * @param length The length of the generated string.
 * @return A random string of the specified length, or an empty string if an error occurs.
 */
std::string generateRandomString(int length) {
#if (defined(__linux__))
    static thread_local std::random_device rd;
    static thread_local std::mt19937 gen(rd());

    std::uniform_int_distribution<> dis('a', 'z');

    std::stringstream ss;
    for (int i = 0; i < length; ++i) {
        ss << static_cast<char>(dis(gen));
    }
    return ss.str();
#else
    HCRYPTPROV hCryptProv;
    if (!CryptAcquireContext(&hCryptProv, NULL, NULL, PROV_RSA_FULL, CRYPT_VERIFYCONTEXT)) {
        std::cerr << "CryptAcquireContext failed, error code: " << GetLastError() << std::endl;
        return "";
    }

    std::stringstream ss;
    BYTE buffer;
    for (int i = 0; i < length; ++i) {
        if (!CryptGenRandom(hCryptProv, sizeof(BYTE), &buffer)) {
            std::cerr << "CryptGenRandom failed, error code: " << GetLastError() << std::endl;
            CryptReleaseContext(hCryptProv, 0);
            return "";
        }
        ss << static_cast<int>(buffer);
    }

    CryptReleaseContext(hCryptProv, 0);
    return ss.str();
#endif
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
void split_and_parallel_align(std::vector<std::string> &data, std::vector<std::string> &name,
                              std::vector<std::vector<std::pair<int_t, int_t>>> &chain) {
    // Print status header
    if (global_args.verbose) {
        std::cout << "#                Parallel Aligning...                       #" << std::endl;
        print_table_divider();
    }

    std::string output = "";
    Timer timer;

    // basic sizes
    uint_t chain_num = 0;
    if (!chain.empty())
        chain_num = chain[0].size();
    uint_t seq_num = data.size();

    // Expand chains first (same as before)
    std::vector<std::vector<std::string>> chain_string(chain_num); // chain_num * seq_num
    std::vector<ExpandChainParams> params(chain_num);
    for (uint_t i = 0; i < chain_num; ++i) {
        params[i].data = &data;
        params[i].chain = &chain;
        params[i].chain_index = i;
        params[i].result_store = chain_string.begin() + i;
    }

    if (global_args.min_seq_coverage == 1) {
        ThreadPool pool(global_args.thread);
        for (uint_t i = 0; i < chain_num; ++i) {
            pool.add_task([&, i]() { expand_chain(&params[i]); });
        }
        pool.shutdown();
    } else {
        for (uint_t i = 0; i < chain_num; ++i) {
            expand_chain(&params[i]);
        }
    }
    params.clear();

    // SW expand timing
    double SW_time = timer.elapsed_time();
    {
        std::stringstream s;
        s << std::fixed << std::setprecision(2) << SW_time;
        if (global_args.verbose) {
            output = "SW expand time: " + s.str() + " seconds.";
            print_table_line(output);
        }
    }

    timer.reset();

    // Calculate parallel alignment ranges and perform parallel alignment for each range
    std::vector<std::vector<std::pair<int_t, int_t>>> parallel_align_range = get_parallel_align_range(data, chain);
    uint_t parallel_num = parallel_align_range.size();

    // We'll produce for each parallel block an MSA built by WFA (center-star)
    std::vector<std::vector<std::string>> parallel_string(parallel_num, std::vector<std::string>(seq_num));

    // Thread pool to run WFA-based MSA for each parallel block
    {
        ThreadPool pool(global_args.thread);

        for (uint_t i = 0; i < parallel_num; ++i) {
            // capture necessary items by value for thread-safety
            auto range = parallel_align_range[i];
            pool.add_task([i, range, &data, &parallel_string]() {
                try {
                    uint_t seq_num_local = data.size();
                    // Build block fragments for each sequence (preserve ordering)
                    std::vector<std::string> block_seqs(seq_num_local);
                    for (uint_t s = 0; s < seq_num_local; ++s) {
                        if (s < range.size()) {
                            auto [start, len] = range[s];
                            if (len > 0 && start >= 0) {
                                const std::string &full = data[s];
                                size_t st = static_cast<size_t>(std::max<int_t>(0, (int_t)start));
                                size_t take = std::min<size_t>((size_t)len, (st < full.size() ? full.size() - st : 0));
                                if (take > 0 && st < full.size())
                                    block_seqs[s] = full.substr(st, take);
                                else
                                    block_seqs[s] = "";
                            } else {
                                block_seqs[s] = "";
                            }
                        } else {
                            // if the range vector is shorter than seq_num, fill with empty
                            block_seqs[s] = "";
                        }
                    }

                    // Call WFA-based MSA (center-star).
                    std::vector<std::string> msa_block = wfa_msa_center_star(block_seqs);

                    // Normalize result: ensure size == seq_num_local
                    if (msa_block.size() != seq_num_local) {
                        // If returned smaller, create vector with gaps and copy what exists
                        size_t final_len = 0;
                        for (const auto &r : msa_block)
                            final_len = std::max(final_len, r.size());
                        if (final_len == 0 && !msa_block.empty())
                            final_len = msa_block[0].size();
                        if (final_len == 0)
                            final_len = 1; // fallback safety

                        std::vector<std::string> normalized(seq_num_local, std::string(final_len, '-'));
                        for (size_t k = 0; k < msa_block.size() && k < seq_num_local; ++k)
                            normalized[k] = msa_block[k];

                        msa_block.swap(normalized);
                    } else {
                        // ensure all rows have same length; if not, pad with gaps
                        size_t L = 0;
                        for (auto &r : msa_block)
                            L = std::max(L, r.size());
                        for (auto &r : msa_block)
                            if (r.size() < L)
                                r.append(L - r.size(), '-');
                    }

                    // store move-constructed result
                    parallel_string[i] = std::move(msa_block);
                } catch (const std::exception &e) {
                    // on error, produce an all-gap block so pipeline can continue
                    size_t safe_len = 1;
                    std::vector<std::string> fallback(data.size(), std::string(safe_len, '-'));
                    parallel_string[i] = std::move(fallback);
                }
            });
        }

        pool.shutdown();
    }
    // FOR SERIAL VERSION (debugging)
    // for (uint_t i = 0; i < parallel_num; ++i) {
    //     auto range = parallel_align_range[i];

    //     try {
    //         uint_t seq_num_local = data.size();
    //         // Build block fragments for each sequence (preserve ordering)
    //         std::vector<std::string> block_seqs(seq_num_local);
    //         for (uint_t s = 0; s < seq_num_local; ++s) {
    //             if (s < range.size()) {
    //                 auto [start, len] = range[s];
    //                 if (len > 0 && start >= 0) {
    //                     const std::string &full = data[s];
    //                     size_t st = static_cast<size_t>(std::max<int_t>(0, (int_t)start));
    //                     size_t take = std::min<size_t>((size_t)len, (st < full.size() ? full.size() - st : 0));
    //                     if (take > 0 && st < full.size())
    //                         block_seqs[s] = full.substr(st, take);
    //                     else
    //                         block_seqs[s] = "";
    //                 } else {
    //                     block_seqs[s] = "";
    //                 }
    //             } else {
    //                 block_seqs[s] = "";
    //             }
    //         }

    //         // Call WFA-based MSA (center-star).
    //         std::vector<std::string> msa_block = wfa_msa_center_star(block_seqs);

    //         // Normalize result: ensure size == seq_num_local
    //         if (msa_block.size() != seq_num_local) {
    //             // If returned smaller, create vector with gaps and copy what exists
    //             size_t final_len = 0;
    //             for (const auto &r : msa_block)
    //                 final_len = std::max(final_len, r.size());
    //             if (final_len == 0 && !msa_block.empty())
    //                 final_len = msa_block[0].size();
    //             if (final_len == 0)
    //                 final_len = 1; // fallback safety

    //             std::vector<std::string> normalized(seq_num_local, std::string(final_len, '-'));
    //             for (size_t k = 0; k < msa_block.size() && k < seq_num_local; ++k)
    //                 normalized[k] = msa_block[k];

    //             msa_block.swap(normalized);
    //         } else {
    //             // ensure all rows have same length; if not, pad with gaps
    //             size_t L = 0;
    //             for (auto &r : msa_block)
    //                 L = std::max(L, r.size());
    //             for (auto &r : msa_block)
    //                 if (r.size() < L)
    //                     r.append(L - r.size(), '-');
    //         }

    //         // store move-constructed result
    //         parallel_string[i] = std::move(msa_block);

    //     } catch (const std::exception &e) {
    //         // on error, produce an all-gap block so pipeline can continue
    //         size_t safe_len = 1;
    //         std::vector<std::string> fallback(data.size(), std::string(safe_len, '-'));
    //         parallel_string[i] = std::move(fallback);
    //     }
    // }

    // Parallel align timing
    double parallel_align_time = timer.elapsed_time();
    {
        std::stringstream s;
        s << std::fixed << std::setprecision(2) << parallel_align_time;
        if (global_args.verbose) {
            output = "Parallel align time: " + s.str() + " seconds.";
            print_table_line(output);
        }
    }

    timer.reset();

    // Concatenate ranges: same as before
    std::vector<std::vector<std::pair<int_t, int_t>>> concat_range = concat_chain_and_parallel_range(chain, parallel_align_range);

    // Normalize each parallel_string[i] to seq_num (filling missing with gaps)
    for (uint_t i = 0; i < parallel_num; ++i) {
        if (parallel_string[i].size() != seq_num) {
            // determine max length in this block and resize/pad
            size_t L = 0;
            for (auto &s : parallel_string[i])
                L = std::max(L, s.size());
            if (L == 0)
                L = 1; // safety
            parallel_string[i].resize(seq_num, std::string(L, '-'));
            for (auto &s : parallel_string[i])
                if (s.size() < L)
                    s.append(L - s.size(), '-');
        }
    }

    // Concatenate chain strings (from expand_chain) and parallel strings (from WFA)
    std::vector<std::vector<std::string>> concat_string = concat_chain_and_parallel(chain_string, parallel_string);
    std::vector<uint_t> fragment_len = get_first_nonzero_lengths(concat_string);

    // seq2profile and concat alignment (same as before)
    seq2profile(concat_string, data, concat_range, fragment_len);
    double seq2profile_time = timer.elapsed_time();

    concat_alignment(concat_string, name);

    // final timings and prints
    {
        std::stringstream s;
        s << std::fixed << std::setprecision(2) << seq2profile_time;
        if (global_args.verbose) {
            output = "Seq-profile time: " + s.str() + " seconds.";
            print_table_line(output);
            print_table_divider();
        }
    }
    return;
}

// parse SAM-style CIGAR (e.g. "10M2I4M")
std::vector<std::pair<int, char>> parse_sam_cigar(const std::string &cigar) {
    std::vector<std::pair<int, char>> ops;
    int num = 0;
    for (unsigned char c : cigar) {
        if (std::isdigit(c)) {
            num = num * 10 + (c - '0');
        } else {
            ops.emplace_back(num, (char)c);
            num = 0;
        }
    }
    return ops;
}

// apply parsed SAM CIGAR to produce two aligned strings (ref, qry)
std::pair<std::string, std::string> apply_cigar_to_seqs(const std::string &cigar, const std::string &ref_seq, const std::string &qry_seq) {

    auto ops = parse_sam_cigar(cigar);
    std::string a_ref;
    a_ref.reserve(ref_seq.size() + qry_seq.size());
    std::string a_qry;
    a_qry.reserve(ref_seq.size() + qry_seq.size());
    size_t i_ref = 0, i_qry = 0;

    for (auto &p : ops) {
        int cnt = p.first;
        char op = p.second;
        if (op == 'M' || op == '=' || op == 'X') {
            for (int k = 0; k < cnt; ++k) {
                if (i_ref >= ref_seq.size() || i_qry >= qry_seq.size())
                    break;
                a_ref.push_back(ref_seq[i_ref++]);
                a_qry.push_back(qry_seq[i_qry++]);
            }
        } else if (op == 'I') {
            // insertion wrt reference => gaps in ref
            for (int k = 0; k < cnt; ++k) {
                if (i_qry >= qry_seq.size())
                    break;
                a_ref.push_back('-');
                a_qry.push_back(qry_seq[i_qry++]);
            }
        } else if (op == 'D' || op == 'N') {
            // deletion wrt reference => gaps in qry
            for (int k = 0; k < cnt; ++k) {
                if (i_ref >= ref_seq.size())
                    break;
                a_ref.push_back(ref_seq[i_ref++]);
                a_qry.push_back('-');
            }
        } else {
            // ignore/defensive: treat as match
            for (int k = 0; k < cnt; ++k) {
                if (i_ref >= ref_seq.size() || i_qry >= qry_seq.size())
                    break;
                a_ref.push_back(ref_seq[i_ref++]);
                a_qry.push_back(qry_seq[i_qry++]);
            }
        }
    }
    return {a_ref, a_qry};
}

// wrapper: run pairwise WFA and return SAM CIGAR string
std::string wfa_pairwise_cigar(const std::string &pattern, const std::string &text, int mismatch, int gap_open, int gap_extend) {
    using namespace wfa; // WFAlignerGapAffine in namespace wfa
    // create an aligner with chosen penalties and default memory mode
    WFAlignerGapAffine aligner(mismatch, gap_open, gap_extend, WFAligner::Alignment, WFAligner::MemoryMed);
    // align end-to-end (pattern = reference, text = query)
    aligner.alignEnd2End(pattern.c_str(), pattern.size(), text.c_str(), text.size());
    // get SAM CIGAR (string). pass false to avoid verbose header
    return aligner.getCIGAR(false);
}

// center-star progressive MSA using WFA pairwise alignments
std::vector<std::string> wfa_msa_center_star(const std::vector<std::string> &sequences) {
    if (sequences.empty())
        return {};
    size_t n = sequences.size();
    if (n == 1)
        return sequences;

    // choose center = longest sequence (simple heuristic)
    size_t center_idx = 0;
    for (size_t i = 1; i < n; ++i)
        if (sequences[i].size() > sequences[center_idx].size())
            center_idx = i;

    // msa_rows holds current aligned rows (may contain gaps)
    std::vector<std::string> msa_rows(n, "");
    msa_rows[center_idx] = sequences[center_idx];

    // align each other sequence to the (progressive) center
    for (size_t idx = 0; idx < n; ++idx) {
        if (idx == center_idx)
            continue;

        // get ungapped version of current center (remove '-' from msa_rows[center_idx])
        std::string ungapped_center;
        ungapped_center.reserve(msa_rows[center_idx].size());
        for (char c : msa_rows[center_idx])
            if (c != '-')
                ungapped_center.push_back(c);

        // pairwise WFA: ungapped_center vs sequences[idx]
        std::string cigar = wfa_pairwise_cigar(ungapped_center, sequences[idx]);
        auto aligned_pair = apply_cigar_to_seqs(cigar, ungapped_center, sequences[idx]);
        std::string new_master = aligned_pair.first;   // may contain '-'
        std::string aligned_qry = aligned_pair.second; // aligned version of sequences[idx]

        // Now we must integrate new_master into the existing msa_rows (propagate gaps)
        // We'll iterate new_master and update every existing row by inserting gaps where new_master has '-'
        std::string &old_master = msa_rows[center_idx]; // existing (may have gaps)
        size_t p = 0;                                   // position in the (mutating) old_master

        // Reconstrói todas as linhas alinhadas de acordo com new_master
        std::vector<std::string> new_rows(n);

        for (size_t i = 0; i < new_master.size(); ++i) {
            char c = new_master[i];
            if (c == '-') {
                // insere gap em todas as linhas
                for (size_t r = 0; r < n; ++r) {
                    new_rows[r].push_back('-');
                }
            } else {
                // avança p até achar posição não-gap em old_master
                while (p < old_master.size() && old_master[p] == '-') {
                    ++p;
                }

                // se p está além do tamanho do old_master, pad com gap
                if (p >= old_master.size()) {
                    for (size_t r = 0; r < n; ++r) {
                        new_rows[r].push_back('-');
                    }
                } else {
                    // copia caractere real da linha antiga
                    for (size_t r = 0; r < n; ++r) {
                        if (p < msa_rows[r].size())
                            new_rows[r].push_back(msa_rows[r][p]);
                        else
                            new_rows[r].push_back('-'); // se a linha acabou, completa com gap
                    }
                }
                ++p;
            }
        }

        // troca as linhas antigas pelas novas
        msa_rows.swap(new_rows);

        // after insertion the center length should equal new_master.size()
        // set center row to new_master (it already has '-' positions consistent with inserted gaps)
        msa_rows[center_idx] = new_master;

        // set aligned query row
        msa_rows[idx] = aligned_qry;

        // pad any existing rows to match new length (defensive)
        size_t L = msa_rows[center_idx].size();
        for (size_t r = 0; r < n; ++r) {
            if (!msa_rows[r].empty() && msa_rows[r].size() < L)
                msa_rows[r].append(L - msa_rows[r].size(), '-');
        }
    }

    // For any rows still empty (shouldn't happen), put original seq and pad
    size_t final_len = msa_rows[center_idx].size();
    std::vector<std::string> out(n);
    for (size_t i = 0; i < n; ++i) {
        if (!msa_rows[i].empty())
            out[i] = msa_rows[i];
        else {
            out[i] = sequences[i];
            if (out[i].size() < final_len)
                out[i].append(final_len - out[i].size(), '-');
        }
    }
    return out;
}

/**
 * @brief Get a vector of integers that are not in the selected_cols vector and have a maximum value of n.
 * @param n The maximum value of the integers in the resulting vector.
 * @param selected_cols A vector of integers that are already selected.
 * @return A vector of integers that are not in selected_cols and have a maximum value of n.
 */
std::vector<int_t> get_remaining_cols(int_t n, const std::vector<int_t> selected_cols) {
    // Initialize a boolean vector to indicate whether an integer is selected.
    std::vector<bool> is_selected(n, false);
    // Mark the integers in the selected_cols vector as selected.
    for (int_t i : selected_cols) {
        if (i < n) {
            is_selected[i] = true;
        }
    }
    // Get the integers that are not selected and store them in a new vector.
    std::vector<int_t> remaining;
    for (int i = 0; i < n; i++) {
        if (!is_selected[i]) {
            remaining.push_back(i);
        }
    }
    return remaining;
}

/**
 * @brief Selects columns from a sequence of split points to enable multi thread.
 * @param split_points_on_sequence A vector of vectors of pairs, where each pair represents the start and mem length
 * @return A vector of indices of the selected columns.
 */
std::vector<int_t> select_columns(std::vector<std::vector<std::pair<int_t, int_t>>> split_points_on_sequence) {
    // Get the number of columns and rows in the split points sequence.
    uint_t col_num = split_points_on_sequence[0].size();
    uint_t row_num = split_points_on_sequence.size();
    // Create a vector to keep track of whether each column needs to be changed.
    std::vector<bool> col_need_change(col_num, false);
    // Create a vector to store the indices of the selected columns.
    std::vector<int_t> selected_cols;
    int_t count = col_num;
    // Create a vector to store the number of effect columns for each column.
    std::vector<std::pair<uint_t, uint_t>> effect_col_num(col_num);

    for (uint_t i = 0; i < col_num; i++) {
        effect_col_num[i].first = 0;
        effect_col_num[i].second = i;
        for (uint_t j = 0; j < row_num; j++) {
            if (split_points_on_sequence[j][i].first == -1) {
                col_need_change[i] = true;
                count--;
                break;
            }
        }
    }

    // Compute the number of effect columns for each column.
    for (uint_t i = 0; i < col_num; i++) {
        if (col_need_change[i]) {
            continue;
        }
        bool has_left = (i == 0 || col_need_change[i - 1] == false);
        bool has_right = (i == col_num - 1 || col_need_change[i + 1] == false);

        if (has_left && has_right) {
            col_need_change[i] = true;
            count--;
            continue;
        }
        effect_col_num[i].first += 1;

        if (has_left && i < col_num - 1) {
            effect_col_num[i + 1].first += 1;
        }

        if (has_right && i > 0) {
            effect_col_num[i - 1].first += 1;
        }
    }
    // Sort the columns based on the number of effect columns.
    std::sort(effect_col_num.begin(), effect_col_num.end(),
              [](const std::pair<uint_t, uint_t> &a, const std::pair<uint_t, uint_t> &b) { return a.first > b.first; });

    for (uint_t i = 0; i < col_num; i++) {
        if (count <= 0) {
            break;
        }
        uint_t effect_num = effect_col_num[i].first;
        uint_t col_index = effect_col_num[i].second;
        if (effect_num <= 0) {
            std::cerr << "some bugs occur in select column." << std::endl;
            exit(-1);
        }
        if (col_need_change[col_index] == false) {
            col_need_change[col_index] = true;
            selected_cols.push_back(col_index);
            count--;
        }
        if (effect_num > 1) {
            if ((col_index == 1 || (col_index >= 2 && col_need_change[col_index - 2]) == true) &&
                (col_index >= 1 && col_need_change[col_index - 1] == false)) {
                col_need_change[col_index - 1] = true;
                count--;
            }

            if ((col_index == col_num - 2 || (col_index + 2 < col_num && col_need_change[col_index + 2]) == true) &&
                (col_index + 1 < col_num && col_need_change[col_index + 1] == false)) {
                col_need_change[col_index + 1] = true;
                count--;
            }
        }
    }
    return selected_cols;
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
    // Cast the input parameters to the correct struct type
    ExpandChainParams *ptr = static_cast<ExpandChainParams *>(arg);
    // Get data, chain, and chain_index from the input parameters
    const std::vector<std::string> data = *(ptr->data);
    std::vector<std::vector<std::pair<int_t, int_t>>> chain = *(ptr->chain);
    const uint_t chain_index = ptr->chain_index;
    // std::cout << "in" << chain_index << '\n';
    // Get the number of sequences in the data vector and the number of chains in the current chain
    uint_t seq_num = data.size();
    uint_t chain_num = chain[0].size();

    // Declares a default Aligner
    StripedSmithWaterman::Aligner aligner;
    // Declares a default filter
    StripedSmithWaterman::Filter filter;
    // Declares an alignment that stores the result
    StripedSmithWaterman::Alignment alignment;

    uint_t query_length = 0;
    std::string query = "";
    std::vector<std::string> aligned_fragment(seq_num);
    // Find the query sequence and its length in the current chain
    for (uint_t i = 0; i < seq_num; i++) {
        if (chain[i][chain_index].first != -1) {
            query_length = chain[i][chain_index].second;
            query = data[i].substr(chain[i][chain_index].first, query_length);
            break;
        }
    }

    for (uint_t i = 0; i < seq_num; i++) {
        int_t begin_pos = chain[i][chain_index].first;
        // If the begin position is -1, the current subsequence is unaligned
        if (begin_pos == -1) {
            uint_t tmp_index = chain_index;
            // Find the beginning and end positions of the unaligned subsequence
            uint_t ref_begin_pos = 0;
            uint_t ref_end_pos = 0;
            int_t maskLen = query_length / 2;
            maskLen = maskLen < 15 ? 15 : maskLen;

            for (; tmp_index > 0 && chain[i][tmp_index - 1].first == -1; --tmp_index)
                ;

            ref_begin_pos = tmp_index <= 0 ? 0 : chain[i][tmp_index - 1].first + chain[i][tmp_index - 1].second;

            tmp_index = chain_index;
            for (; tmp_index < chain_num - 1 && chain[i][tmp_index + 1].first == -1; ++tmp_index)
                ;

            ref_end_pos = tmp_index >= chain_num - 1 ? data[i].length() - 1 : chain[i][tmp_index + 1].first;

            std::string ref = data[i].substr(ref_begin_pos, ref_end_pos - ref_begin_pos);

            // Get the reference subsequence and align it with the query subsequence
            aligner.Align(query.c_str(), ref.c_str(), ref.size(), filter, &alignment, maskLen);

            std::pair<int_t, int_t> p = store_sw_alignment(alignment, ref, query, aligned_fragment, i);

            if (p.first != -1) {
                p.first += ref_begin_pos;
                (*(ptr->chain))[i][chain_index] = p;
            }

        } else {
            aligned_fragment[i] = query;
        }
    }
    *(ptr->result_store) = aligned_fragment;

    return NULL;
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
    // Extract cigar string from the alignment
    std::vector<unsigned int> cigar = alignment.cigar;
    // Extract the start and end positions of the alignment on the reference sequence
    int_t ref_begin = alignment.ref_begin;
    int_t ref_end = alignment.ref_end;

    uint_t query_begin = 0;

    std::string aligned_result = "";
    // If the alignment failed, return (-1,-1)
    if (ref_begin <= -1) {
        res_store[seq_index] = aligned_result;
        return std::make_pair(-1, -1);
    }

    int_t S_count = 0;
    int_t total_length = alignment.query_end;
    int_t new_ref_begin = ref_begin;
    uint_t new_ref_end = ref_end;
    if (cigar_int_to_op(cigar[0]) == 'S') {
        S_count += cigar_int_to_len(cigar[0]);
        if (new_ref_begin - cigar_int_to_len(cigar[0]) < 0) {
            new_ref_begin = 0;
        } else {
            new_ref_begin = new_ref_begin - cigar_int_to_len(cigar[0]);
        }
    }
    if (cigar_int_to_op(cigar[cigar.size() - 1]) == 'S') {
        S_count += cigar_int_to_len(cigar[cigar.size() - 1]);
        total_length += cigar_int_to_len(cigar[cigar.size() - 1]);
        if (new_ref_end + cigar_int_to_len(cigar[cigar.size() - 1]) > (uint_t)ref.length() - 1) {
            new_ref_end = ref.length() - 1;
        } else {
            new_ref_end = new_ref_end + cigar_int_to_len(cigar[cigar.size() - 1]);
        }
    }

    if (S_count > ceil(0.8 * total_length)) {
        res_store[seq_index] = aligned_result;
        return std::make_pair(-1, -1);
    }
    uint_t p_ref = 0;
    uint_t p_query = 0;

    for (uint_t i = 0; i < cigar.size(); i++) {
        char op = cigar_int_to_op(cigar[i]);
        int_t len = cigar_int_to_len(cigar[i]);

        switch (op) {
        case 'S': {
            // Handle soft clipping at the beginning and end of the alignment
            if (i == 0) {
                int_t tmp_len = len;
                if (ref_begin <= len) {
                    while (tmp_len > ref_begin) {
                        aligned_result += "-";
                        tmp_len--;
                    }
                }
                for (int_t j = ref_begin - tmp_len; j < ref_begin; j++) {
                    aligned_result += ref[j];
                }
            } else {
                int_t tmp_len = len;
                for (uint_t j = ref_end + 1; j < ref.length() && tmp_len > 0; j++) {
                    aligned_result += ref[j];
                    tmp_len--;
                }
                while (tmp_len > 0) {
                    aligned_result += "-";
                    tmp_len--;
                }
            }
            p_query += len;
            break;
        }
        case 'M':
        case 'X':
        case '=': {
            // Handle match, mismatch, and substitution operations
            for (int_t j = 0; j < len; j++) {
                aligned_result += ref[ref_begin + p_ref + j];
            }
            p_ref += len;
            p_query += len;
            break;
        }
        case 'I': {
            // Handle insertion operations
            for (int_t j = 0; j < len; j++) {
                aligned_result += '-';
            }
            p_query += len;
            break;
        }
        case 'D': {
            // Handle deletion operations
            std::string gaps = "";
            for (int_t j = 0; j < len; j++) {
                gaps += '-';
            }
            for (int_t j = 0; j < len; j++) {
                aligned_result += ref[ref_begin + p_ref + j];
            }
            p_ref += len;

            query.insert(query_begin + p_query, gaps);
            for (uint_t j = 0; j < seq_index; j++) {
                if (res_store[j].length() != 0) {
                    res_store[j].insert(query_begin + p_query, gaps);
                }
            }
            p_query += len;
            break;
        }

        default:
            break;
        }
    }

    res_store[seq_index] = aligned_result;

    std::pair<int, int> p(new_ref_begin, new_ref_end - new_ref_begin + 1);
    return p;
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
void concat_alignment(std::vector<std::vector<std::string>> &concat_string, std::vector<std::string> &name) {
    std::string output_path = global_args.output_path;
    std::vector<std::string> concated_data(name.size(), "");
    // Concatenate the sequences
    for (uint_t i = 0; i < name.size(); i++) {
        for (uint_t j = 0; j < concat_string.size(); j++) {
            concated_data[i] += concat_string[j][i];
        }
    }
    // Write the concatenated sequences to the output file
    std::ofstream output_file;
    output_file.open(output_path);
    if (!output_file.is_open()) {
        std::cerr << "Error opening output file " << output_path << std::endl;
        exit(1);
    }

    for (uint_t i = 0; i < concated_data.size(); i++) {
        std::stringstream ss;
        ss << ">" << name[i] << "\n" << concated_data[i] << "\n";
        output_file << ss.str();
    }
    output_file.close();
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
