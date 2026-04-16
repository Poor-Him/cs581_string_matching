#!/usr/bin/env python3
import argparse
import os
import time


FASTA_EXTENSIONS = {".fa", ".fasta", ".fna", ".ffn", ".faa", ".frn"}


def build_prefix(pattern_idx):
    m = len(pattern_idx)
    pi = [0] * m
    k = 0

    for q in range(1, m):
        while k > 0 and pattern_idx[k] != pattern_idx[q]:
            k = pi[k - 1]
        if pattern_idx[k] == pattern_idx[q]:
            k += 1
        pi[q] = k

    return pi


def build_dfa(pattern, alphabet):
    char_to_idx = {ch: i for i, ch in enumerate(alphabet)}
    pattern_idx = [char_to_idx[ch] for ch in pattern]

    m = len(pattern)
    s = len(alphabet)

    if m == 0:
        return [[0] * s], char_to_idx

    pi = build_prefix(pattern_idx)

    dfa = [[0] * s for _ in range(m + 1)]
    dfa[0][pattern_idx[0]] = 1

    for q in range(1, m + 1):
        fallback = pi[q - 1]
        dfa[q][:] = dfa[fallback][:]
        if q < m:
            dfa[q][pattern_idx[q]] = q + 1

    return dfa, char_to_idx


def detect_file_type(path, requested_type):
    if requested_type != "auto":
        return requested_type

    ext = os.path.splitext(path)[1].lower()
    if ext in FASTA_EXTENSIONS:
        return "fasta"

    with open(path, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            stripped = line.lstrip()
            if not stripped:
                continue
            return "fasta" if stripped.startswith(">") else "text"

    return "text"


def iter_file_characters(path, file_type, ignore_case=False):
    """
    Yields characters from the input file in a streaming fashion.

    fasta mode:
      - skips header lines starting with '>'
      - ignores whitespace/newlines inside sequence lines

    text mode:
      - yields every character, including spaces/newlines
    """
    with open(path, "r", encoding="utf-8", errors="ignore") as f:
        if file_type == "fasta":
            for line in f:
                if line.startswith(">"):
                    continue

                if ignore_case:
                    line = line.upper()

                for ch in line:
                    if ch.isspace():
                        continue
                    yield ch

        else:  # text mode
            for chunk in f:
                if ignore_case:
                    chunk = chunk.upper()

                for ch in chunk:
                    yield ch


def count_matches_in_file(path, file_type, pattern_len, dfa, char_to_idx, ignore_case=False):
    state = 0
    count = 0
    seq_len = 0

    for ch in iter_file_characters(path, file_type, ignore_case=ignore_case):
        seq_len += 1
        idx = char_to_idx.get(ch, -1)
        state = 0 if idx == -1 else dfa[state][idx]

        if state == pattern_len:
            count += 1

    return count, seq_len


def main():
    parser = argparse.ArgumentParser(
        description="DFA string matcher for FASTA/FNA genomes and normal text files"
    )
    parser.add_argument(
        "--input-file", "--sequence-file",
        dest="input_file",
        required=True,
        help="Path to input file (.txt, .fasta, .fna, etc.)"
    )
    parser.add_argument(
        "--pattern-file",
        default=None,
        help="Read pattern from a file instead of --pattern"
    )
    parser.add_argument(
        "--pattern",
        default=None, 
        help="Pattern string to search for"
    )
    parser.add_argument(
        "--file-type",
        choices=["auto", "fasta", "text"],
        default="auto",
        help="Input type: auto-detect, fasta, or text"
    )
    parser.add_argument(
        "--alphabet",
        default=None,
        help=(
            "Alphabet to use for DFA construction. "
            "If omitted: FASTA uses ACGTN, text uses unique chars from the pattern."
        )
    )
    parser.add_argument(
        "--ignore-case",
        action="store_true",
        help="Convert both pattern and file contents to uppercase before matching"
    )
    parser.add_argument(
        "--repeat",
        type=int,
        default=1,
        help="Repeat matching and keep the best time"
    )

    args = parser.parse_args()

    file_type = detect_file_type(args.input_file, args.file_type)

    if args.pattern_file:
        with open(args.pattern_file, encoding="utf-8") as f:
            pattern = f.read().strip()
    else:
        pattern = args.pattern.replace("\n", "\\n")
    
    if args.ignore_case:
        pattern = pattern.upper()

    if not pattern:
        raise ValueError("Pattern must not be empty.")

    if args.alphabet is not None:
        alphabet = list(dict.fromkeys(args.alphabet.upper() if args.ignore_case else args.alphabet))
    else:
        if file_type == "fasta":
            alphabet = list("ACGTN")
        else:
            alphabet = list(dict.fromkeys(pattern))

    alphabet_set = set(alphabet)

    missing = set(pattern) - alphabet_set
    if missing:
        raise ValueError(
            f"Pattern contains characters not in alphabet: {sorted(missing)}"
        )

    t0 = time.perf_counter()
    dfa, char_to_idx = build_dfa(pattern, alphabet)
    t1 = time.perf_counter()

    best_time = float("inf")
    best_count = None
    best_len = None

    total_time = 0

    for _ in range(args.repeat):
        start = time.perf_counter()
        count, seq_len = count_matches_in_file(
            args.input_file,
            file_type,
            len(pattern),
            dfa,
            char_to_idx,
            ignore_case=args.ignore_case
        )
        end = time.perf_counter()

        elapsed = end - start
        total_time += elapsed
        if elapsed < best_time:
            best_time = elapsed
            best_count = count
            best_len = seq_len

    preprocess_time = t1 - t0

    print(f"Input file:         {args.input_file}")
    print(f"Detected type:      {file_type}")
    print(f"Pattern length (m): {len(pattern)}")
    print(f"Input length (n):   {best_len}")
    print(f"Alphabet size:      {len(alphabet)}")
    print(f"Preprocessing time: {preprocess_time:.8f} s")
    print(f"Matching time:      {best_time:.8f} s")
    print(f"Total time:         {preprocess_time + best_time:.8f} s")
    print(f"Average time:       {(total_time/args.repeat):.8f} s")
    print(f"Match count:        {best_count}")


if __name__ == "__main__":
    main()