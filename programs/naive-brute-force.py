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

def naive(text, pattern):
    #print("Number of words:", len(text.split()))
    #print("Number of unique words:", len(set(text.split())))
    n = len(text)
    m = len(pattern)
    validshifts = 0
    startMTime = time.perf_counter()
    for s in range(0, n-m+1):
        patternMatch = True
        for i in range(m):
            if(text[s+i] != pattern[i]):
                patternMatch = False
                break
        if(patternMatch):
            validshifts += 1
            #print(f"Valid shift {s}")
    #print(f"Valid shifts: {validshifts}")
    endMTime = time.perf_counter()
    return validshifts, endMTime-startMTime

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

if args.ignore_case:
    pattern = pattern.upper()
if not pattern:
    raise ValueError("Pattern must not be empty.")

print(f"Input file:         {args.input_file}")
print(f"Detected type:      {file_type}")
print(f"Pattern length (m): {len(pattern)}")
text = loadFile(args.input_file, file_type)
print(f"Input length (n):   {len(text)}")
best_time = float('inf')
for _ in range(args.repeat):
    match_count, match_time = naive(text, pattern)
    best_time = min(best_time, match_time)
match_time = best_time
#print(f"Time taken: {(end-start):.9f}\n")

print(f"Alphabet size:      {len(alphabet)}")
print(f"Preprocessing time: {0:.8f} s")
print(f"Matching time:      {(match_time):.8f} s")
print(f"Total time:         {(match_time):.8f} s")
print(f"Match count:        {match_count}")

'''
for item in [
    ["text1.txt", "AW"],
    ["dna1000.txt", "gaa"],
    ["newsarticle.txt", "help"]
]:
    start = time.perf_counter()
    naive(item[0], item[1])
    print(f"Time taken: {(time.perf_counter()-start):.9f}\n")
'''