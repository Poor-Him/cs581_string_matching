#!/usr/bin/env python3
import time
import os
import argparse

FASTA_EXTENSIONS = {".fa", ".fasta", ".fna", ".ffn", ".faa", ".frn"}

def loadFile(filename, filetype):
    text = ""
    # Get this directory
    base_dir = os.path.dirname(os.path.abspath(__file__))
    # Navigate up one level, then into Texts
    texts_dir = os.path.join(base_dir, "..", "Texts")
    if(filetype == "text"):
        with open(os.path.join(texts_dir, filename), encoding="utf-8") as f:
            for line in f:
                text += line.strip()
                text += "\\n"
    else:
        with open(os.path.join(texts_dir, filename), encoding="utf-8") as f:
            for line in f:
                if line.startswith(">"):
                    continue
                text += line.strip()
    #print(text)
    return text

def rabin_karp(text, pattern, d, q):
    n = len(text)
    m = len(pattern)
    h = pow(d, m-1, q)
    p = 0
    ts = 0
    validShifts = 0
    startPTime = time.perf_counter()
    for i in range(m):
        p = (d*p+ ord(pattern[i])) % q
        ts = (d*ts+ ord(text[i])) % q
    endPTime = time.perf_counter()
    startMTime = time.perf_counter()
    for s in range(0, n-m+1):
        if p == ts:
            patternMatch = True
            for i in range(m):
                if(text[s+i] != pattern[i]):
                    patternMatch = False
                    break
            if(patternMatch):
                validShifts += 1
        if(s < n-m):
            ts = (d*(ts - ord(text[s]) * h) + ord(text[s+m])) % q
    #print(f"Valid shifts: {validShifts}")
    endMTime = time.perf_counter()
    return validShifts, endPTime-startPTime, endMTime-startMTime


parser = argparse.ArgumentParser(
    description="Naive Brute-force string matching algorithm"
)
parser.add_argument(
    "--input-file", "--sequence-file",
    dest="input_file",
    required=True,
    help="Path to input file (.txt or .fna)"
)
parser.add_argument(
    "--pattern",
    required=True,
    help="Pattern string to search for"
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
    "--file-type",
    choices=["auto", "fasta", "text"],
    default="auto",
    help="Input type: auto-detect, fasta, or text"
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

file_type = args.file_type
ext = os.path.splitext(args.input_file)[1].lower()
if(ext == ".fna"):
    file_type = "fasta"
elif(ext == ".txt"):
    file_type = "text"
else:
    file_type = "unknown"

pattern = args.pattern.replace("\n", "")
if args.alphabet is not None:
    alphabet = list(dict.fromkeys(args.alphabet.upper() if args.ignore_case else args.alphabet))
else:
    if file_type == "fasta":
        alphabet = list("ACGTN")
    else:
        alphabet = list(dict.fromkeys(pattern))
alphabet_set = set(alphabet)


d = 256
q = 1000003


if args.ignore_case:
    pattern = pattern.upper()
if not pattern:
    raise ValueError("Pattern must not be empty.")

print(f"Input file:         {args.input_file}")
print(f"Detected type:      {file_type}")
print(f"Pattern length (m): {len(pattern)}")
text = loadFile(args.input_file, file_type)
print(f"Input length (n):   {len(text)}")
bestptime = float('inf')
bestmtime = float('inf')
for _ in range(args.repeat):
    match_count, ptime, mtime = rabin_karp(text, pattern, d, q)
    if(bestptime + bestmtime > ptime+mtime):
        bestptime = ptime
        bestmtime = mtime
#print(f"Time taken: {(end-start):.9f}\n")

print(f"Alphabet size:      {len(alphabet)}")
print(f"Preprocessing time: {bestptime:.8f} s")
print(f"Matching time:      {(bestmtime):.8f} s")
print(f"Total time:         {(bestptime + bestmtime):.8f} s")
print(f"Match count:        {match_count}")