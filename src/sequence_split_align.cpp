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
    // === CONCATENATE CHAIN + PARALLEL ===
    auto concat_blocks = concat_chain_and_parallel(chain_string, spoa_parallel_string);

    // === FINAL OUTPUT ===
    concat_alignment_from_blocks(concat_blocks, name);

    if (global_args.verbose) {
        print_table_line(std::format("Seq-profile time: {:.2f} seconds.", timer.elapsed_time()));
        print_table_divider();
    }
}

/**
 * @brief Concatena blocos alinhados em uma única saída FASTA
 *
 * Esta função é a **ponte final** entre o alinhamento por blocos
 * e a saída. Ela toma os blocos já alinhados (MEM + INTERVALO)
 * e os concatena em sequências completas.
 *
 * @param blocks Vector de blocos alinhados (cada bloco = MSA de um intervalo/MEM)
 * @param names Nomes das sequências
 * @return void (escreve em arquivo)
 *
 * Nota: Remplaces a chamada antiga a:
 *   concat_alignment(concatstring, name);
 *
 * Mas agora com fluxo simplificado:
 *   concat_alignment_from_blocks(concat_blocks, name);
 */
void concat_alignment_from_blocks(const std::vector<std::vector<std::string>> &blocks, const std::vector<std::string> &names) {

    // ========================================================================
    // Validação de Entrada
    // ========================================================================

    if (blocks.empty()) {
        if (global_args.verbose) {
            print_table_line("Warning: Empty blocks for concatenation");
        }
        return;
    }

    size_t seq_num = names.size();

    if (seq_num == 0) {
        throw std::runtime_error("Error: No sequence names provided");
    }

    // Validar que cada bloco tem o número correto de sequências
    for (size_t b = 0; b < blocks.size(); ++b) {
        if (blocks[b].size() != seq_num) {
            throw std::runtime_error(std::format("Error: Block {} has {} sequences, expected {}", b, blocks[b].size(), seq_num));
        }
    }

    // ========================================================================
    // Concatenação de Blocos
    // ========================================================================

    // Pré-alocar strings para evitar reallocations repetidas
    std::vector<std::string> final_alignment(seq_num);

    // Estimar tamanho total para reserva eficiente
    size_t total_length = 0;
    for (size_t b = 0; b < blocks.size(); ++b) {
        for (size_t s = 0; s < seq_num; ++s) {
            total_length += blocks[b][s].size();
        }
    }

    // Reservar memória para evitar reallocations
    for (size_t s = 0; s < seq_num; ++s) {
        final_alignment[s].reserve(total_length / seq_num + 100);
    }

    // Concatenar blocos em ordem
    for (size_t b = 0; b < blocks.size(); ++b) {
        const auto &block = blocks[b];

        for (size_t s = 0; s < seq_num; ++s) {
            // Verificar se bloco[s] não está vazio antes de concatenar
            if (!block[s].empty()) {
                final_alignment[s] += block[s];
            }
        }
    }

    // ========================================================================
    // Validação de Alinhamento (após concatenação)
    // ========================================================================

    // Verificar que todos têm o mesmo comprimento (propriedade MSA)
    size_t first_len = final_alignment[0].size();

    for (size_t s = 1; s < seq_num; ++s) {
        if (final_alignment[s].size() != first_len) {
            throw std::runtime_error(std::format("Error: After concatenation, sequence {} has length {} "
                                                 "but sequence 0 has length {}",
                                                 s, final_alignment[s].size(), first_len));
        }
    }

    if (global_args.verbose) {
        print_table_line(std::format("Final alignment: {} sequences × {} bp", seq_num, first_len));
    }

    // ========================================================================
    // Escrever Saída FASTA
    // ========================================================================

    std::string output_path = global_args.output_path;
    std::ofstream output_file(output_path);

    if (!output_file.is_open()) {
        throw std::runtime_error(std::format("Error: Cannot open output file '{}'", output_path));
    }

    // Escrever cada sequência em formato FASTA
    for (size_t s = 0; s < seq_num; ++s) {
        // Cabeçalho FASTA
        output_file << ">" << names[s] << "\n";

        // Sequência alinhada (com possíveis line breaks)
        const std::string &seq = final_alignment[s];
        const size_t LINE_WIDTH = 80; // Quebra linha a cada 80 caracteres

        for (size_t i = 0; i < seq.size(); i += LINE_WIDTH) {
            size_t chunk_size = std::min(LINE_WIDTH, seq.size() - i);
            output_file << seq.substr(i, chunk_size) << "\n";
        }
    }

    output_file.close();

    if (global_args.verbose) {
        print_table_line(std::format("Output written to: {}", output_path));
    }
}

/**
 * @brief Extrai fragmentos baseado em ranges
 * @param data Sequências originais
 * @param range Pares (start, length) para cada sequência
 * @param fragments Output: fragmentos extraídos
 * @return Comprimento máximo entre os fragmentos
 */
size_t extract_fragments(const std::vector<std::string> &data, const std::vector<std::pair<int, int>> &range,
                         std::vector<std::string> &fragments) {
    size_t max_len = 0;
    fragments.resize(data.size());

    for (uint_t s = 0; s < data.size(); ++s) {
        // ✓ Validação: verificar bounds
        if (s >= range.size()) {
            fragments[s] = "";
            continue;
        }

        auto [start, len] = range[s];

        // ✓ Validação: verificar valores
        if (start < 0 || len <= 0 || start >= (int)data[s].size()) {
            fragments[s] = "";
            continue;
        }

        // ✓ Validação: não ultrapassar tamanho
        if (start + len > (int)data[s].size()) {
            len = data[s].size() - start;
        }

        fragments[s] = data[s].substr(start, len);
        max_len = std::max(max_len, fragments[s].size());
    }

    return max_len;
}

/**
 * @brief Verifica se todos os fragmentos são iguais
 */
bool all_fragments_equal(const std::vector<std::string> &fragments) {
    if (fragments.empty() || fragments[0].empty())
        return false;
    const auto &first = fragments[0];
    return std::all_of(fragments.begin(), fragments.end(), [&first](const auto &s) { return s == first; });
}

/**
 * @brief Extrai sub-fragmento de um fragmento
 */
std::string extract_sub_fragment(const std::string &fragment, size_t pos, size_t end) {
    // ✓ Validação: verificar bounds
    if (fragment.empty() || pos >= fragment.size()) {
        return "";
    }

    // ✓ Garantir que pos < end
    if (pos >= end) {
        return "";
    }

    size_t real_start = std::min(pos, fragment.size());
    size_t real_end = std::min(end, fragment.size());

    if (real_start < real_end) {
        return fragment.substr(real_start, real_end - real_start);
    }

    return "";
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
        }
    }

    // Run the SPOA alignment on the extracted fragments
    // and store the result at the designated location
    *(params->result_store) = spoa_align(fragments);

    return nullptr;
}

/**
 * @brief Merge simples: concatena sub-blocos alinhados cortando overlap
 *
 * Estratégia: NÃO re-alinha! Apenas corta o overlap do fim do bloco anterior
 * e cola o bloco atual.
 *
 * @param sub_results Alinhamentos dos sub-blocos
 * @param original Fragmentos originais (para referência)
 * @param overlap Tamanho do overlap
 * @return Alinhamento completo do bloco
 */
std::vector<std::string> merge_subdivisions_simple(const std::vector<std::vector<std::string>> &sub_results,
                                                   const std::vector<std::string> &original, size_t overlap) {

    size_t seq_num = original.size();
    std::vector<std::string> merged(seq_num, "");

    // ✓ Validação: se sub_results vazio, retorna original
    if (sub_results.empty()) {
        return original;
    }

    // ✓ Validação: primeira entrada pode estar vazia
    if (sub_results[0].size() != seq_num) {
        // Sub_results corrompido, retorna original
        return original;
    }

    // First sub-block: copy completely
    merged = sub_results[0];

    // Subsequent sub-blocks
    for (size_t sub_idx = 1; sub_idx < sub_results.size(); ++sub_idx) {
        // ✓ Validação: verificar tamanho
        if (sub_results[sub_idx].size() != seq_num) {
            continue; // Skip corrupted sub-result
        }

        const auto &sub = sub_results[sub_idx];

        // For each sequence
        for (size_t seq = 0; seq < seq_num; ++seq) {
            // ✓ Validação: verificar se referências são válidas
            if (seq >= merged.size() || seq >= sub.size()) {
                continue;
            }

            const std::string &prev = merged[seq];
            const std::string &curr = sub[seq];

            // ✓ Validação: curr pode ser vazio
            if (curr.empty()) {
                continue;
            }

            // Calculate effective overlap
            size_t overlap_len = std::min(overlap, prev.size());

            // Remove overlap from end of merged
            if (overlap_len > 0 && prev.size() >= overlap_len) {
                merged[seq].erase(prev.size() - overlap_len);
            }

            // Concatenate current sub-block
            merged[seq] += curr;
        }
    }

    return merged;
}

/**
 * @brief Versão simplificada: SPOA apenas, sem duplo alinhamento
 *
 * Estratégia:
 * 1. Blocos iguais → cópia direta
 * 2. Blocos pequenos (<10 kbp) → SPOA direto
 * 3. Blocos grandes (≥10 kbp) → subdivide + align + merge_simples
 *
 * @param data Sequências originais
 * @param parallel_align_range Ranges dos intervalos
 * @param pool Thread pool
 * @return Alinhamentos SPOA prontos (sem reprocessamento)
 */
std::vector<std::vector<std::string>> preprocess_parallel_blocks(const std::vector<std::string> &data,
                                                                 const std::vector<std::vector<std::pair<int, int>>> &parallel_align_range,
                                                                 ThreadPool &pool) {

    const size_t MAX_BLOCK_SIZE = 10000;
    const size_t OVERLAP_SIZE = 400;

    uint_t parallel_num = parallel_align_range.size();
    uint_t seq_num = data.size();

    std::vector<std::vector<std::string>> result(parallel_num, std::vector<std::string>(seq_num));

    std::atomic<uint_t> count_exact(0);
    std::atomic<uint_t> count_spoa_direct(0);
    std::atomic<uint_t> count_spoa_subdivided(0);

    const auto *data_ptr = &data;
    const auto *range_ptr = &parallel_align_range;

    // ========================================================================
    // PASS 1: Schedule ALL blocks as tasks
    // ========================================================================

    for (uint_t block_idx = 0; block_idx < parallel_num; ++block_idx) {
        pool.add_task([block_idx, data_ptr, range_ptr, &result, &count_exact, &count_spoa_direct, &count_spoa_subdivided]() {
            const auto &data = *data_ptr;
            const auto &parallel_align_range_ref = *range_ptr;

            // ✓ Validação: verificar índice
            if (block_idx >= parallel_align_range_ref.size()) {
                return;
            }

            const auto &range = parallel_align_range_ref[block_idx];

            // Extract fragments for this block
            std::vector<std::string> fragments(data.size());
            size_t max_len = extract_fragments(data, range, fragments);

            // === CASE 1: All fragments equal? -> Copy directly ===
            if (all_fragments_equal(fragments)) {
                std::fill(result[block_idx].begin(), result[block_idx].end(), fragments[0]);
                count_exact.fetch_add(1, std::memory_order_relaxed);
                return;
            }

            // === CASE 2: Small block (<10kbp)? -> Direct SPOA alignment ===
            if (max_len <= MAX_BLOCK_SIZE) {
                result[block_idx] = run_spoa_local(fragments);
                count_spoa_direct.fetch_add(1, std::memory_order_relaxed);
                return;
            }

            // === CASE 3: Large block (≥10kbp) -> Subdivide + Align + Merge ===
            count_spoa_subdivided.fetch_add(1, std::memory_order_relaxed);

            size_t stride = MAX_BLOCK_SIZE - OVERLAP_SIZE;
            size_t num_subblocks = (max_len + stride - 1) / stride;

            // ✓ Validação: se num_subblocks é 0, não processar
            if (num_subblocks == 0) {
                result[block_idx] = fragments;
                return;
            }

            std::vector<std::vector<std::string>> sub_results(num_subblocks);

            size_t sub_idx = 0;
            for (size_t pos = 0; pos < max_len; pos += stride) {
                size_t end = std::min(pos + MAX_BLOCK_SIZE + OVERLAP_SIZE, max_len);

                std::vector<std::string> sub_frags(data.size());
                for (uint_t s = 0; s < data.size(); ++s) {
                    sub_frags[s] = extract_sub_fragment(fragments[s], pos, end);
                }

                // ✓ Try-catch para SPOA que pode falhar
                try {
                    sub_results[sub_idx] = run_spoa_local(sub_frags);
                } catch (...) {
                    // Se SPOA falha, usa original
                    sub_results[sub_idx] = sub_frags;
                }

                if (end >= max_len)
                    break;

                ++sub_idx;
            }

            // ✓ Merge com sub_results que pode estar parcialmente vazio
            result[block_idx] = merge_subdivisions_simple(sub_results, fragments, OVERLAP_SIZE);
        });
    }

    // ========================================================================
    // PASS 2: Wait for ALL blocks to complete
    // ========================================================================
    pool.wait_for_tasks();

    if (global_args.verbose) {
        print_table_line("Blocos iguais (cópia): " + std::to_string(count_exact.load(std::memory_order_relaxed)));
        print_table_line("Blocos SPOA direto: " + std::to_string(count_spoa_direct.load(std::memory_order_relaxed)));
        print_table_line("Blocos SPOA subdivido: " + std::to_string(count_spoa_subdivided.load(std::memory_order_relaxed)));
    }

    return result;
}

/**
 * @brief Runs SPOA alignment on a set of sequences (helper for subdivisions)
 * @param seqs Input sequences (may contain empty strings)
 * @return Aligned sequences with gaps
 */
std::vector<std::string> run_spoa_local(const std::vector<std::string> &seqs) {
    if (seqs.empty())
        return {};

    try {
        // Use SPOA parameters (match=5, mismatch=-4, gap=-8)
        auto aligner = spoa::AlignmentEngine::Create(spoa::AlignmentType::kNW, 5, -4, -8);
        spoa::Graph graph;

        std::vector<size_t> non_empty_indices;
        non_empty_indices.reserve(seqs.size());

        // Add only non-empty sequences
        for (size_t i = 0; i < seqs.size(); ++i) {
            if (seqs[i].empty())
                continue;

            try {
                auto alignment = aligner->Align(seqs[i], graph);
                graph.AddAlignment(alignment, seqs[i]);
                non_empty_indices.push_back(i);
            } catch (...) {
                // Se alguma sequência causar problema, skip
                continue;
            }
        }

        // If all empty, return empty strings
        if (non_empty_indices.empty()) {
            return std::vector<std::string>(seqs.size(), "");
        }

        // Generate MSA
        auto msa_compact = graph.GenerateMultipleSequenceAlignment();
        size_t aln_len = msa_compact.empty() ? 0 : msa_compact[0].size();

        // ✓ Validação: se MSA vazio, retorna original
        if (aln_len == 0) {
            return seqs;
        }

        // Rebuild to original size: empty sequences become all gaps
        std::vector<std::string> msa(seqs.size(), std::string(aln_len, '-'));
        for (size_t k = 0; k < non_empty_indices.size(); ++k) {
            if (k < msa_compact.size()) {
                msa[non_empty_indices[k]] = std::move(msa_compact[k]);
            }
        }

        return msa;
    } catch (...) {
        // Se SPOA falha completamente, retorna original
        return seqs;
    }
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
 * @brief Concatena MEMs com intervalos alinhados (SEM refinement)
 *
 * Estratégia:
 * - Intercala: MEM[0] + INTERVALO[0] + MEM[1] + INTERVALO[1] + ...
 * - NÃO chama refinement() (perde informação!)
 * - Cada bloco é retornado como está
 * - Concatenação final feita após esta função
 *
 * @param chain_string MEMs alinhados (do expand_chain)
 * @param parallel_string Intervalos alinhados (do preprocess_parallel_blocks_v2)
 * @return Blocos intercalados MEM+INTERVALO
 */
std::vector<std::vector<std::string>> concat_chain_and_parallel_v2(std::vector<std::vector<std::string>> &chain_string,
                                                                   std::vector<std::vector<std::string>> &parallel_string) {

    size_t seq_num = chain_string.empty() ? 0 : chain_string[0].size();
    if (seq_num == 0 && !parallel_string.empty()) {
        seq_num = parallel_string[0].size();
    }

    size_t chain_num = chain_string.size();
    size_t parallel_num = parallel_string.size();

    std::vector<std::vector<std::string>> result;
    result.reserve(chain_num + parallel_num);

    // Intercalar: INTERVALO[0] + MEM[0] + INTERVALO[1] + MEM[1] + ...
    // (ordem pode variar conforme seu get_parallel_align_range)

    for (size_t i = 0; i < parallel_num; ++i) {
        result.push_back(parallel_string[i]); // INTERVALO
        if (i < chain_num) {
            result.push_back(chain_string[i]); // MEM
        }
    }

    // Blocos MEM restantes (se houver)
    for (size_t i = parallel_num; i < chain_num; ++i) {
        result.push_back(chain_string[i]);
    }

    // ✗ NÃO chama refinement()!
    // ✗ NÃO chama concat_alignment() aqui!
    // Apenas retorna os blocos intercalados

    return result;
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
    // Get maximum size for chains (number of chain blocks)
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
            std::pair<int_t, int_t> chain_entry = {-1, -1};
            if (j < chain[i].size()) {
                chain_entry = chain[i][j];
            } else {
                // If this sequence has fewer chain entries, treat missing entries as gaps
                chain_entry = {-1, -1};
            }
            int_t begin_pos = chain_entry.first;

            if (begin_pos == -1) {
                tmp_range.push_back({-1, -1});
                last_pos = -1;
            } else {
                if (last_pos == -1) {
                    // previous was missing, so current also mark as missing for the parallel fragment
                    tmp_range.push_back({-1, -1});
                } else {
                    // normal non-overlapping region: from last_pos length (begin_pos - last_pos)
                    tmp_range.push_back({last_pos, begin_pos - last_pos});
                }
                // update last_pos safely (ensure not negative)
                if (chain_entry.second >= 0)
                    last_pos = begin_pos + chain_entry.second;
                else
                    last_pos = begin_pos; // fallback, shouldn't happen
            }
        }

        // Handle the end of the sequence after the last chain
        if (last_pos == -1) {
            tmp_range.push_back({-1, -1});
        } else {
            int_t remaining = static_cast<int_t>(data[i].length()) - last_pos;
            if (remaining < 0)
                remaining = 0;
            tmp_range.push_back({last_pos, remaining});
        }
    }

    // Allocate transpose result vector
    std::vector<std::vector<std::pair<int_t, int_t>>> transpose_res;
    transpose_res.reserve(chain_num + 1);

    // Transpose parallel_align_range into transpose_res
    for (uint_t j = 0; j <= chain_num; j++) {
        std::vector<std::pair<int_t, int_t>> transposed_row(seq_num);
        for (uint_t i = 0; i < seq_num; i++) {
            // safe read: if j is out of range for this row, set {-1,-1}
            if (j < parallel_align_range[i].size())
                transposed_row[i] = parallel_align_range[i][j];
            else
                transposed_row[i] = {-1, -1};
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
 * @brief Concatenate two sets of sequence data (chain and parallel) into a single set of concatenated data.
 * @param chain_string A vector of vectors containing the chain sequence data.
 * @param parallel_string A vector of vectors containing the parallel sequence data.
 * @return std::vector<std::vectorstd::string> A vector of vectors containing the concatenated sequence data.
 */
std::vector<std::vector<std::string>> concat_chain_and_parallel(std::vector<std::vector<std::string>> &chain_string,
                                                                std::vector<std::vector<std::string>> &parallel_string) {
    // Handle empty inputs
    if (parallel_string.empty() && chain_string.empty())
        return {};
    uint_t seq_num = 0;
    if (!parallel_string.empty())
        seq_num = parallel_string[0].size();
    else if (!chain_string.empty())
        seq_num = chain_string[0].size();

    uint_t chain_num = chain_string.size();
    uint_t parallel_num = parallel_string.size();
    std::vector<std::vector<std::string>> concated_data(chain_num + std::max<uint_t>(parallel_num, 1));
    uint_t count = 0;

    for (uint_t i = 0; i < parallel_num; i++) {
        if (parallel_string[i].size() != seq_num) {
            // defensive: resize if needed
            for (auto &s : parallel_string[i]) {
                if (s.size() != seq_num) {
                    // ignore, will handle later
                }
            }
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
    // if there are leftover chain_string entries and no parallel entries
    if (parallel_num == 0) {
        for (uint_t i = 0; i < chain_num; ++i) {
            for (uint_t j = 0; j < seq_num; ++j) {
                concated_data[i].push_back(chain_string[i][j]);
            }
        }
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
