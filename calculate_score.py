#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Unified, fast, dependency-free MSA scoring:
- SP (Sum-of-Pairs), Avg SP, Scaled SP
- Entropy (per column summed)
- Percentage of Non-Gaps
- Percentage of Totally Conserved Columns
- Average pairwise identity (computed column-wise)
Parallelized with multiprocessing (fork-friendly).
Author: ChatGPT Assistant - Nov/2025
Modeled by: Thiago Luiz Parolin
"""

from __future__ import annotations
import argparse
import math
import os
from multiprocessing import Pool, cpu_count
from typing import List, Tuple

# -------------------------
# FAST FASTA read + preprocess
# -------------------------
def read_fasta(path: str) -> Tuple[List[str], List[str]]:
    ids = []
    seqs = []
    with open(path, 'r') as f:
        curr_id = None
        buf = []
        for line in f:
            line = line.rstrip('\n')
            if not line:
                continue
            if line[0] == '>':
                if curr_id is not None:
                    seqs.append(''.join(buf))
                curr_id = line[1:].split()[0]
                ids.append(curr_id)
                buf = []
            else:
                buf.append(line.strip())
        if curr_id is not None:
            seqs.append(''.join(buf))
    # uppercase and normalize U->T
    seqs = [s.upper().replace('U', 'T') for s in seqs]
    return ids, seqs

def preprocess_replace_invalid(seqs: List[str]) -> List[str]:
    # replace any char not ACTGN- with N
    out = []
    valid = set(b'ACGTN-')
    for s in seqs:
        b = bytearray(s.encode('ascii', 'ignore'))
        for i in range(len(b)):
            c = b[i]
            if c not in valid:
                b[i] = ord('N')
        out.append(b.decode('ascii'))
    return out

# -------------------------
# Worker globals (set in initializer)
# -------------------------
SEQ_BYTES = None           # list of bytes objects
NSEQ = 0
LSEQ = 0
LOOKUP = None              # array[256] -> small int index
MATCH_S = 2.0
MISMATCH_S = -1.0
GAP1_S = -2.0
GAP2_S = 0.0

# mapping indices:
# 0:A,1:C,2:G,3:T,4:N,5:-
# names for clarity
IDX_A, IDX_C, IDX_G, IDX_T, IDX_N, IDX_DASH = range(6)

# initializer for pool workers
def _init_worker(seq_bytes, matchS, mismatchS, gap1S, gap2S):
    global SEQ_BYTES, NSEQ, LSEQ, LOOKUP, MATCH_S, MISMATCH_S, GAP1_S, GAP2_S
    SEQ_BYTES = seq_bytes
    NSEQ = len(SEQ_BYTES)
    LSEQ = len(SEQ_BYTES[0])
    MATCH_S = matchS
    MISMATCH_S = mismatchS
    GAP1_S = gap1S
    GAP2_S = gap2S
    # build lookup table
    tbl = [IDX_N] * 256
    tbl[ord('A')] = IDX_A
    tbl[ord('C')] = IDX_C
    tbl[ord('G')] = IDX_G
    tbl[ord('T')] = IDX_T
    tbl[ord('N')] = IDX_N
    tbl[ord('-')] = IDX_DASH
    # also accept lowercase if any (shouldn't after preprocess)
    tbl[ord('a')] = IDX_A
    tbl[ord('c')] = IDX_C
    tbl[ord('g')] = IDX_G
    tbl[ord('t')] = IDX_T
    tbl[ord('n')] = IDX_N
    tbl[ord('-')] = IDX_DASH
    LOOKUP = tbl

# -------------------------
# core: process block of columns [start, end)
# returns aggregated tuple of metrics
# -------------------------
def _process_block(args):
    start, end = args
    # localize globals
    seqs = SEQ_BYTES
    lookup = LOOKUP
    nseq = NSEQ

    total_sp = 0.0
    total_entropy = 0.0
    non_gap_bases = 0     # total count of non-gap tokens across all columns
    totally_conserved_cols = 0
    pairwise_matches_sum = 0  # sum over columns of number of matching pairs
    pairwise_comparable_sum = 0  # sum over columns of number of comparable pairs (n_non_gap choose 2)
    cols_processed = 0

    for j in range(start, end):
        # count vector for this column
        c0=c1=c2=c3=c4=c5 = 0  # A C G T N dash
        for s in seqs:
            v = lookup[s[j]]
            if v == IDX_A:
                c0 += 1
            elif v == IDX_C:
                c1 += 1
            elif v == IDX_G:
                c2 += 1
            elif v == IDX_T:
                c3 += 1
            elif v == IDX_N:
                c4 += 1
            else:
                c5 += 1

        # SP-like counts as in original script
        A = c0; C = c1; G = c2; T = c3; N = c4; dash = c5
        # match pairs: sum comb(count[x],2) for A,C,G,T
        match_pairs = (A*(A-1) + C*(C-1) + G*(G-1) + T*(T-1)) // 2
        # mismatch as original (keeps same semantics)
        mismatch_pairs = ((A + C) * (G + T)) + A*C + G*T
        gap1_pairs = (A + C + G + T + N) * dash
        gap2_pairs = (dash * (dash - 1) // 2) + ((N * (N - 1)) // 2) + (A + C + G + T) * N

        col_sp = (match_pairs * MATCH_S) + (mismatch_pairs * MISMATCH_S) + (gap1_pairs * GAP1_S) + (gap2_pairs * GAP2_S)
        total_sp += col_sp

        # pairwise identity aggregates:
        # number of matching pairs among non-gaps:
        match_pairs_non_gaps = (A*(A-1) + C*(C-1) + G*(G-1) + T*(T-1) + N*(N-1)) // 2
        # comparable pairs = comb(non_gaps, 2)
        non_gaps = A + C + G + T + N
        comparable_pairs = (non_gaps * (non_gaps - 1)) // 2
        pairwise_matches_sum += match_pairs_non_gaps
        pairwise_comparable_sum += comparable_pairs

        # entropy over symbols A,C,G,T,N (exclude gaps) — if no non-gaps, entropy=0
        if non_gaps > 0:
            # compute probabilities
            # use safe log2
            ent = 0.0
            for cnt in (A, C, G, T, N):
                if cnt > 0:
                    p = cnt / non_gaps
                    ent -= p * math.log2(p)
            total_entropy += ent
        else:
            total_entropy += 0.0

        # non-gap base count for "percentage non-gaps"
        non_gap_bases += non_gaps

        # totally conserved: all non-gap same and no gaps and no N
        if dash == 0 and N == 0:
            # exactly one non-zero among A,C,G,T and equals nseq
            if ((A == nseq) or (C == nseq) or (G == nseq) or (T == nseq)):
                totally_conserved_cols += 1

        cols_processed += 1

    return (total_sp, total_entropy, non_gap_bases, totally_conserved_cols,
            pairwise_matches_sum, pairwise_comparable_sum, cols_processed)

# -------------------------
# formatting (Brazilian)
# -------------------------
def fmt_int_br(x: int) -> str:
    s = f"{x:,}"
    # return s.replace(',', '.')
    return s

def fmt_float_br(x: float, ndigits: int = 4) -> str:
    # round then replace decimal dot with comma and thousand separator with dot
    fmt = f"{{:,.{ndigits}f}}".format(x)
    # return fmt.replace(',', 'X').replace('.', ',').replace('X', '.')
    return fmt

# -------------------------
# main driver
# -------------------------
def main():
    parser = argparse.ArgumentParser(description="Fast MSA scoring (SP, entropy, coverage, conserved cols).")
    parser.add_argument('--input', '-i', type=str, required=True, help='Aligned FASTA file')
    parser.add_argument('--match', type=float, default=2.0, help='match score')
    parser.add_argument('--mismatch', type=float, default=-1.0, help='mismatch score')
    parser.add_argument('--gap1', type=float, default=-2.0, help='gap-base score')
    parser.add_argument('--gap2', type=float, default=0.0, help='gap-gap score')
    parser.add_argument('--procs', '-p', type=int, default=max(1, cpu_count()-1), help='Number of worker processes')
    parser.add_argument('--digits', type=int, default=4, help='Decimal digits for floats in output')
    args = parser.parse_args()

    ids, seqs = read_fasta(args.input)
    if len(seqs) == 0:
        print("No sequences found.")
        return
    if len(seqs) == 1:
        print("Only one sequence found; metrics require at least 2 sequences.")
        return

    # preprocess sequences (replace invalid chars with N)
    seqs = preprocess_replace_invalid(seqs)

    # sanity: ensure all same length
    L = len(seqs[0])
    for s in seqs:
        if len(s) != L:
            raise ValueError("All sequences must have the same length (aligned).")

    # convert to list of bytes for fast indexing
    seq_bytes = [s.encode('ascii') for s in seqs]
    nseq = len(seq_bytes)

    # prepare pool and split columns into chunks
    procs = max(1, min(args.procs, cpu_count()))
    # create chunk boundaries (equal sized)
    chunk_size = max(1, (L + procs - 1) // procs)
    blocks = []
    for i in range(0, L, chunk_size):
        blocks.append((i, min(L, i + chunk_size)))

    # initialize pool (fork-friendly: prepare seq_bytes before Pool)
    pool = Pool(processes=procs, initializer=_init_worker,
                initargs=(seq_bytes, args.match, args.mismatch, args.gap1, args.gap2))

    results = pool.map(_process_block, blocks)
    pool.close()
    pool.join()

    # aggregate
    total_sp = total_entropy = non_gap_bases = totally_conserved_cols = 0
    pairwise_matches_sum = pairwise_comparable_sum = cols_processed = 0
    for r in results:
        (sp, ent, ngb, tcc, pms, pcs, cp) = r
        total_sp += sp
        total_entropy += ent
        non_gap_bases += ngb
        totally_conserved_cols += tcc
        pairwise_matches_sum += pms
        pairwise_comparable_sum += pcs
        cols_processed += cp

    # final metrics
    total_pairs = (nseq * (nseq - 1)) // 2
    avg_sp = total_sp / total_pairs if total_pairs > 0 else 0.0
    scaled_sp = avg_sp / L if L > 0 else 0.0
    mean_entropy = total_entropy / cols_processed if cols_processed > 0 else 0.0
    pct_non_gaps = (non_gap_bases / (nseq * L)) * 100.0
    pct_totally_conserved = (totally_conserved_cols / L) * 100.0

    # average pairwise identity (column-wise)
    avg_pairwise_identity = (pairwise_matches_sum / pairwise_comparable_sum) if pairwise_comparable_sum > 0 else 0.0

    # output (Brazilian format)
    print("=== MSA scoring results ===")
    print(f"Sequences: {fmt_int_br(nseq)}")
    print(f"Columns : {fmt_int_br(L)}")
    print()
    print(f"SP score (sum-of-pairs) : {fmt_float_br(total_sp, args.digits)}")
    print(f"Avg SP score (per pair) : {fmt_float_br(avg_sp, args.digits)}")
    print(f"Scaled SP score         : {fmt_float_br(scaled_sp, args.digits)}")
    print()
    print(f"Entropy (mean per col)  : {fmt_float_br(mean_entropy, args.digits)} bits")
    print(f"Percentage non-gaps     : {fmt_float_br(pct_non_gaps, 2)} %")
    print(f"Totally conserved cols  : {fmt_float_br(pct_totally_conserved, 4)} %")
    print(f"Average pairwise identity (col-wise) : {fmt_float_br(avg_pairwise_identity, 6)}")
    print("===========================")

if __name__ == "__main__":
    main()
