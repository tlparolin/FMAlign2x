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

#include "mem_finder.h"

#if defined(M64) && M64 == 1
#define LIBSAIS_OMP libsais64_omp
#define LIBSAIS_PLCP_OMP libsais64_plcp_omp
#else
#define LIBSAIS_OMP libsais_omp
#define LIBSAIS_PLCP_OMP libsais_plcp_omp
#endif

void *find_optimal_chain(void *arg) {
    auto *ptr = static_cast<FindOptimalChainParams *>(arg);
    auto &chains = *ptr->chains;
    const size_t chain_num = chains.size();

    std::vector<double> dp(chain_num);
    std::vector<size_t> prev(chain_num, static_cast<size_t>(-1));
    // Iterate over all "mem" objects and calculate their size and update dynamic programming tables
    for (size_t i = 0; i < chain_num; ++i) {
        dp[i] = static_cast<double>(chains[i].second);
        for (size_t j = i + 1; j < chain_num; ++j) {
            if (chains[i].first + chains[i].second < chains[j].first && dp[i] > dp[j]) {
                dp[j] = dp[i];
                prev[j] = i;
            }
        }
    }
    // Find the index of the last "mem" object in the longest non-conflicting sequence
    double max_size = 0;
    size_t end_index = 0;
    for (size_t i = 0; i < chain_num; ++i) {
        if (dp[i] > max_size) {
            max_size = dp[i];
            end_index = i;
        }
    }
    // Retrieve the indices of all non-conflicting "mem" objects in the longest sequence
    std::vector<bool> selected(chain_num, false);
    while (end_index < chain_num && end_index != static_cast<size_t>(-1)) {
        selected[end_index] = true;
        end_index = prev[end_index];
    }

    std::vector<std::pair<int_t, int_t>> new_chains;
    new_chains.reserve(chain_num);
    for (size_t i = 0; i < chain_num; ++i) {
        if (selected[i])
            new_chains.emplace_back(chains[i]);
        else
            new_chains.emplace_back(-1, -1);
    }

    chains = std::move(new_chains);

    return nullptr;
}

/**
 * @brief DP sequence number times! Filter out overlapping memory regions and generate split points for each sequence.
 * Given a vector of memory regions and the number of sequences, this function removes any
 * overlapping memory regions and generates split points for each sequence based on the non-overlapping regions.
 * @param mems Vector of memory regions.
 * @param sequence_num Number of sequences.
 * @return Vector of split points for each sequence.
 */
std::vector<std::vector<std::pair<int_t, int_t>>> filter_mem_accurate(std::vector<mem> &mems, size_t sequence_num) {
    // delete MEM full of "-"
    std::erase_if(mems, [](const mem &m) { return m.mem_length <= 0; });

    const size_t mem_num = mems.size();

    // Initialize a vector of vectors of pairs of integers to represent the split points for each sequence
    std::vector<std::vector<std::pair<int_t, int_t>>> split_point_on_sequence(sequence_num,
                                                                              std::vector<std::pair<int_t, int_t>>(mem_num, {-1, -1}));

    // Loop through each non-conflicting MEM in the input
    for (size_t i = 0; i < mem_num; ++i) {
        // Get the current MEM and its substring positions
        const auto &m = mems[i];
        // Loop through each substring of the current MEM
        for (const auto &substring : m.substrings) {
            // Create a pair of the substring position and the length of the MEM
            auto &current_point = split_point_on_sequence[substring.sequence_index][i];
            std::pair<int_t, int_t> candidate = {substring.position, m.mem_length};
            // If this split point is already set for this sequence and it is farther from the average position,
            // skip this split point and move to the next one
            if (current_point.first != -1) {
                if (std::abs(candidate.first - m.avg_pos) > std::abs(current_point.first - m.avg_pos)) {
                    continue;
                }
            }
            // Set this split point for this sequence to the current substring position and MEM length
            current_point = candidate;
        }
    }

    std::vector<FindOptimalChainParams> find_params(sequence_num);
    ThreadPool pool(global_args.thread);
    for (size_t i = 0; i < sequence_num; ++i) {
        find_params[i].chains = &split_point_on_sequence[i];
        pool.add_task([&, i]() { find_optimal_chain(&find_params[i]); });
    }
    pool.shutdown(); // Waits for all tasks to finish and finalizes the pool

    // remove column that too much -1
    std::vector<size_t> selected_cols;
    const size_t col_count = split_point_on_sequence[0].size();
    const size_t max_missing = static_cast<size_t>(std::floor(sequence_num * (1.0 - global_args.min_seq_coverage)));

    for (size_t j = 0; j < col_count; ++j) {
        size_t missing = 0;
        for (size_t i = 0; i < sequence_num; ++i) {
            if (split_point_on_sequence[i][j].first == -1) {
                ++missing;
            }
        }
        if (missing <= max_missing) {
            selected_cols.push_back(j);
        }
    }

    std::vector<std::vector<std::pair<int_t, int_t>>> chain(sequence_num);
    for (size_t i = 0; i < sequence_num; ++i) {
        chain[i].reserve(selected_cols.size());
        for (size_t col : selected_cols) {
            chain[i].emplace_back(split_point_on_sequence[i][col]);
        }
    }

    return chain;
}

/**
 * @brief DP only Once!Filter out overlapping memory regions and generate split points for each sequence.
 * Given a vector of memory regions and the number of sequences, this function removes any
 * overlapping memory regions and generates split points for each sequence based on the non-overlapping regions.
 * @param mems Vector of memory regions.
 * @param sequence_num Number of sequences.
 * @return Vector of split points for each sequence.
 */
std::vector<std::vector<std::pair<int_t, int_t>>> filter_mem_fast(std::vector<mem> &mems, uint_t sequence_num) {
    // delete MEM full of "-"
    std::vector<mem>::iterator mem_it = mems.begin();
    while (mem_it != mems.end()) {
        mem tmp_mem = *mem_it;
        if (tmp_mem.mem_length <= 0) {
            mem_it = mems.erase(mem_it);
        } else {
            mem_it++;
        }
    }
    // Initialize dynamic programming tables to keep track of size and previous indices
    uint_t mem_num = mems.size();
    std::vector<double> dp(mem_num, 0);
    std::vector<int_t> prev(mem_num, -1);
    // Iterate over all "mem" objects and calculate their size and update dynamic programming tables
    for (uint_t i = 0; i < mem_num; ++i) {
        double size = mems[i].mem_length * mems[i].substrings.size();
        dp[i] += size;
        for (uint_t j = i + 1; j < mem_num; j++) {
            if (mems[i].avg_pos + mems[i].mem_length < mems[j].avg_pos && dp[i] > dp[j]) {
                dp[j] = dp[i];
                prev[j] = i;
            }
        }
    }
    // Find the index of the last "mem" object in the longest non-conflicting sequence
    double max_size = 0;
    int_t end_index = 0;
    for (uint_t i = 0; i < mem_num; i++) {
        if (dp[i] > max_size) {
            max_size = dp[i];
            end_index = i;
        }
    }
    // Retrieve the indices of all non-conflicting "mem" objects in the longest sequence
    std::vector<int_t> mems_without_conflict;
    while (end_index > 0) {
        mems_without_conflict.push_back(end_index);
        end_index = prev[end_index];
    }
    reverse(mems_without_conflict.begin(), mems_without_conflict.end());
    // Initialize a vector of vectors of pairs of integers to represent the split points for each sequence
    std::vector<std::vector<std::pair<int_t, int_t>>> split_point_on_sequence(
        sequence_num, std::vector<std::pair<int_t, int_t>>(mems_without_conflict.size(), std::make_pair(-1, -1)));

    // Loop through each non-conflicting MEM in the input
    for (uint_t i = 0; i < mems_without_conflict.size(); i++) {
        // Get the current MEM and its substring positions
        mem tmp_mem = mems[mems_without_conflict[i]];
        // Loop through each substring of the current MEM
        for (uint_t j = 0; j < tmp_mem.substrings.size(); j++) {
            // Create a pair of the substring position and the length of the MEM
            std::pair<int_t, int_t> p(tmp_mem.substrings[j].position, tmp_mem.mem_length);
            // If this split point is already set for this sequence and it is farther from the average position,
            // skip this split point and move to the next one
            if (split_point_on_sequence[tmp_mem.substrings[j].sequence_index][i].first != -1) {
                if (abs(p.first - tmp_mem.avg_pos) >
                    abs(split_point_on_sequence[tmp_mem.substrings[j].sequence_index][i].first - tmp_mem.avg_pos)) {
                    continue;
                }
            }
            // Set this split point for this sequence to the current substring position and MEM length
            split_point_on_sequence[tmp_mem.substrings[j].sequence_index][i] = p;
        }
    }

    // Loop through each sequence in the input
    for (uint_t i = 0; i < sequence_num; i++) {
        // Initialize the index of the last split point on this sequence to 0
        int_t last_end_index = 0;
        // Loop through each pair of split points on this sequence that do not conflict
        for (uint_t j = 1; j < mems_without_conflict.size(); j++) {
            // Get the position and length of the current split point
            int_t cur_pos = split_point_on_sequence[i][j].first;
            int_t cur_len = split_point_on_sequence[i][j].second;
            // If the current split point has a negative position, skip it
            if (cur_pos < 0) {
                continue;
            }
            if (cur_pos >= split_point_on_sequence[i][last_end_index].first &&
                (j == split_point_on_sequence[0].size() - 1 || cur_pos <= split_point_on_sequence[i][j + 1].first)) {
                if (cur_pos >= split_point_on_sequence[i][last_end_index].first + split_point_on_sequence[i][last_end_index].second) {
                    last_end_index = j;
                }
                // If the current split point conflicts with the last split point used,
                // choose the split point with the shortest length and mark the other one as invalid
                else {
                    if (split_point_on_sequence[i][last_end_index].second > cur_len) {
                        split_point_on_sequence[i][j].first = -1;
                        split_point_on_sequence[i][j].second = -1;
                    } else {
                        split_point_on_sequence[i][last_end_index].first = -1;
                        split_point_on_sequence[i][last_end_index].second = -1;
                        last_end_index = j;
                    }
                }
            } else {
                split_point_on_sequence[i][j].first = -1;
                split_point_on_sequence[i][j].second = -1;
            }
            // If the current split point is after the last split point that was used,
            // update the index of the last split point used to the current index
        }
    }
    // remove column that too much -1
    std::vector<int_t> selected_cols;
    for (uint_t j = 0; j < split_point_on_sequence[0].size(); j++) {
        int_t count = 0;
        for (uint_t i = 0; i < split_point_on_sequence.size(); i++) {
            if (split_point_on_sequence[i][j].first == -1) {
                count++;
            }
        }
        if (count <= floor(sequence_num * (1 - global_args.min_seq_coverage))) {
            selected_cols.push_back(j);
        }
    }
    std::vector<std::vector<std::pair<int_t, int_t>>> chain(sequence_num);
    for (uint_t i = 0; i < selected_cols.size(); i++) {
        for (uint_t j = 0; j < split_point_on_sequence.size(); j++) {
            chain[j].push_back(split_point_on_sequence[j][selected_cols[i]]);
        }
    }
    return chain;
}

/**
 * @brief Find MEMs in a set of sequences.
 * @param data A vector of strings representing the sequences.
 * @return Vector of split points for each sequence.
 */
std::vector<std::vector<std::pair<int_t, int_t>>> find_mem(const std::vector<std::string> &data) {
    if (global_args.verbose) {
        std::cout << "#                    Finding MEM...                         #" << std::endl;
        print_table_divider();
    }

    Timer timer;
    size_t n = 0;

    unsigned char *concat_data = concat_strings(data, n);

    if (global_args.min_mem_length < 0) {
        int_t l = ceil(pow(n, 1 / (global_args.degree + 2)));
        l = l > 30 ? l : 30;
        l = l < 2000 ? l : 2000;
        global_args.min_mem_length = l;
    }

    if (global_args.verbose) {
        print_table_line("Minimal MEM length is set to " + std::to_string(global_args.min_mem_length));
    }

    if (global_args.filter_mode == "default") {
        global_args.filter_mode = (data.size() < 100) ? "local" : "global";
    }

    if (global_args.verbose) {
        print_table_line("Filter mode is set to " + global_args.filter_mode);
    }

    if (global_args.min_seq_coverage < 0) {
        global_args.min_seq_coverage = (data.size() < 100) ? 1 : 0.7;
    }

    if (global_args.verbose) {
        print_table_line("Minimal sequence coverage is set to " + std::to_string(global_args.min_seq_coverage));
    }

    timer.reset();

    std::vector<int_t> SA(n);
    std::vector<int_t> PLCP(n);

    LIBSAIS_OMP(concat_data, SA.data(), n, 0, NULL, global_args.thread);
    LIBSAIS_PLCP_OMP(concat_data, SA.data(), PLCP.data(), n, global_args.thread);

    double suffix_construction_time = timer.elapsed_time();
    std::stringstream s;
    s << std::fixed << std::setprecision(2) << suffix_construction_time;
    if (global_args.verbose) {
        print_table_line("Suffix construction time: " + s.str() + " seconds");
    }

    timer.reset();
    const int_t min_mem_length = global_args.min_mem_length;
    const int_t min_cross_sequence = std::ceil(global_args.min_seq_coverage * data.size());

    // Build a vector of joined sequence boundaries.
    std::vector<uint_t> joined_sequence_bound;
    joined_sequence_bound.reserve(data.size());
    uint_t total_length = 0;
    for (const auto &seq : data) {
        joined_sequence_bound.push_back(total_length);
        total_length += seq.length() + 1;
    }
    auto intervals = get_lcp_intervals(PLCP.data(), SA.data(), min_mem_length, min_cross_sequence, n);
    const uint_t interval_size = intervals.size();

    // Allocate space for MEMs and conversion parameters.
    std::vector<mem> mems(interval_size);
    // Convert each interval to a MEM in parallel
    IntervalToMemConversionParams *params = new IntervalToMemConversionParams[interval_size];

    ThreadPool pool(global_args.thread);
    for (uint_t i = 0; i < interval_size; i++) {
        params[i].SA = &SA;
        params[i].interval = intervals[i];
        params[i].concat_data = concat_data;
        params[i].result_store = mems.begin() + i;
        params[i].min_mem_length = min_mem_length;
        params[i].joined_sequence_bound = joined_sequence_bound;

        pool.add_task([&, i]() { interval2mem(params + i); });
    }
    pool.shutdown();

    if (mems.empty() && global_args.verbose) {
        print_table_line("Warning: There is no MEMs, please adjust your parameters.");
    }

    // Sort the MEMs and map them to the original data.
    sort_mem(mems, data);

    free(concat_data);
    delete[] params;

    uint_t sequence_num = data.size();
    auto split_point_on_sequence =
        (global_args.filter_mode == "global") ? filter_mem_fast(mems, sequence_num) : filter_mem_accurate(mems, sequence_num);

    global_args.avg_file_size = (n / (split_point_on_sequence[0].size() + 1)) / std::pow(2, 20);

    double mem_time = timer.elapsed_time();
    std::stringstream ss;
    ss << std::fixed << std::setprecision(3) << mem_time;
    if (global_args.verbose) {
        print_table_line("Sequence divide parts (MEMs): " + std::to_string(split_point_on_sequence[0].size() + 1));
        print_table_line("MEM process time: " + ss.str() + " seconds");
        print_table_divider();
    }

    return split_point_on_sequence;
}

/**
 * @brief Concatenates a vector of strings with separator 1 and a terminating 0.
 * @param strings The vector of strings to concatenate.
 * @param n A reference to the total length of the concatenated string.
 * @return A pointer to the concatenated string.
 * @note The returned string must be deleted by the caller.
 */
unsigned char *concat_strings(const std::vector<std::string> &strings, size_t &n) {
    // Calculate total size required (including separators and terminator)
    n = std::accumulate(strings.begin(), strings.end(), std::size_t(1),
                        [](std::size_t sum, const std::string &s) { return sum + s.size() + 1; });

    // Allocate buffer
    auto concat_data = std::make_unique<unsigned char[]>(n);

    // Copy strings to buffer
    std::size_t index = 0;
    for (const auto &s : strings) {
        std::memcpy(concat_data.get() + index, s.data(), s.size());
        index += s.size();
        concat_data[index++] = 1; // Separador
    }

    // Add terminator 0
    concat_data[n - 1] = 0;

    return concat_data.release(); // Transferência de propriedade para o chamador
}

/**
 * @brief an LCP (Longest Common Prefix) array and a threshold value,
 * finds all the LCP intervals where each value is greater than or equal to the threshold value,
 * and at least one value in the interval is equal to the threshold value.
 * @param lcp_array The input LCP array
 * @param threshold The threshold value
 * @param min_cross_sequence the min number of crossed sequence
 * @return  The output vector of pairs representing the LCP intervals
 */
std::vector<std::pair<uint_t, uint_t>> get_lcp_intervals(const int_t *plcp_array, int_t *sa, int_t threshold, int_t min_cross_sequence,
                                                         uint_t n) {
    std::vector<std::pair<uint_t, uint_t>> intervals;
    intervals.reserve(n / (min_cross_sequence > 0 ? min_cross_sequence : 1) + 1);

    if (global_args.verbose) {
        print_table_line("Minimal cross sequence number: " + std::to_string(min_cross_sequence));
    }

    int_t left = 0, right = 0;
    bool found = false;

    while (right < (int_t)n) {
        int_t current_val = plcp_array[sa[right]];
        if (current_val >= threshold) { // Changed from LCP to PLCP[SA]
            if (current_val == threshold) {
                found = true;
            }
            right++;
        } else {
            if (found && right - left + 1 >= min_cross_sequence) {
                intervals.emplace_back(left, right);
            }
            ++right;
            left = right;
            found = false;
        }
    }

    if (found && right - left + 1 >= min_cross_sequence) {
        intervals.emplace_back(left, right);
    }
    return intervals;
}

/**
 *@brief This function converts an LCP interval to a MEM (Maximal Exact Match).
 *@param arg A void pointer to the input parameters.
 *@return void* A void pointer to the result, which is stored in the input parameters structure.
 */
void *interval2mem(void *arg) {
    IntervalToMemConversionParams *ptr = static_cast<IntervalToMemConversionParams *>(arg);

    const int_t *SA = ptr->SA->data();
    const int_t min_mem_length = ptr->min_mem_length;
    const unsigned char *concat_data = ptr->concat_data;
    const std::vector<uint_t> &joined_sequence_bound = ptr->joined_sequence_bound;
    std::pair<uint_t, uint_t> interval = ptr->interval;

    mem result;
    uint_t *mem_index = new uint_t(0);
    result.mem_index = mem_index;
    result.mem_length = 0;
    std::vector<uint_t> mem_position;
    mem_position.reserve(interval.second - interval.first + 1); // Reservar espaço para evitar realocações

    for (uint_t i = interval.first - 1; i < interval.second; i++) {
        sub_string tmp_substring;

        auto it = std::upper_bound(joined_sequence_bound.begin(), joined_sequence_bound.end(), SA[i]);
        uint_t sequence_index = std::distance(joined_sequence_bound.begin(), it) - 1;

        tmp_substring.sequence_index = sequence_index;
        tmp_substring.position = SA[i] - joined_sequence_bound[sequence_index];
        mem_position.push_back(SA[i]);

        tmp_substring.mem_index = mem_index;
        result.substrings.push_back(tmp_substring);
    }

    uint_t offset = 1;
    bool all_char_same = true;

    while (true) {
        if (mem_position[0] < offset)
            break;
        char current_char = concat_data[mem_position[0] - offset];

        all_char_same = true;
        for (uint_t i = 1; i < mem_position.size(); i++) {
            if (mem_position[i] < offset || current_char != concat_data[mem_position[i] - offset]) {
                all_char_same = false;
                break;
            }
        }

        if (!all_char_same)
            break;
        offset++;
    }

    if (offset > 0)
        offset--;

    result.mem_length = min_mem_length + offset;
    for (auto &substring : result.substrings) {
        substring.position -= offset;
    }

    uint_t gap_count = 0;
    uint_t base_position = result.substrings[0].position + joined_sequence_bound[result.substrings[0].sequence_index];
    for (int_t i = 0; i < result.mem_length; i++) {
        if (concat_data[base_position + i] == '-') {
            gap_count++;
        }
    }

    if (gap_count > ceil(0.8 * result.mem_length)) {
        result.mem_length = -1;
    }

    *(ptr->result_store) = result;

    return NULL;
}

// function to compute average position of a mem in sequences
void compute_mem_avg_pos(mem &m) {
    float_t sum_pos = 0;
    for (const auto &s : m.substrings) {
        sum_pos += s.position;
    }
    m.avg_pos = sum_pos / m.substrings.size();
}

/**
 *Sorts the input vector of MEMs by the average position of each MEM's substrings along the sequences.
 *Removes any MEMs that span across multiple sequences.
 *Assigns a unique index to each MEM based on its position in the sorted vector.
 *@param mems The vector of MEMs to be sorted.
 *@param data The vector of sequences used to compute the MEMs.
 */
void sort_mem(std::vector<mem> &mems, const std::vector<std::string> &data) {
    // Remove MEMs inválidos e computa a posição média em uma única passagem
    auto it = std::remove_if(mems.begin(), mems.end(), [&](mem &m) {
        if (m.substrings[0].position + m.mem_length >= data[m.substrings[0].sequence_index].length()) {
            return true;
        }
        compute_mem_avg_pos(m);
        return false;
    });
    mems.erase(it, mems.end());

    // Ordena os MEMs pela posição média
    std::sort(std::execution::par, mems.begin(), mems.end(), [](const mem &m1, const mem &m2) { return m1.avg_pos < m2.avg_pos; });

    // Atribui mem_index com base na posição no vetor ordenado
    for (uint_t i = 0; i < mems.size(); i++) {
        *mems[i].mem_index = i;
    }
}
