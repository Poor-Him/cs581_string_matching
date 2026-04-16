
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

def constructpi(P, pi):
    start_preproc = time.perf_counter()
    length = 0
    m = len(P)
    pi[0] = 0
    q = 1
    while q < m:
        if P[q] == P[length]:
            length += 1
            pi[q] = length
            q += 1
        elif length != 0:
            length = pi[length - 1]
        else:
            pi[q] = 0
            q += 1
    end_preproc = time.perf_counter()
    return pi, end_preproc - start_preproc

def search(P, T):
    start_search = time.perf_counter()
    n = len(T)
    m = len(P)
    pi = [0] * m
    res = []
    pi, preproc_time = constructpi(P, pi)
    i = 0
    j = 0
    while i < n:
        if T[i] == P[j]:
            i += 1
            j += 1
            if j == m:
                res.append(i-j)
                j = pi[j-1]
        elif j != 0:
            j = pi[j-1]
        else:
            i += 1
    end_search = time.perf_counter()
    search_time = end_search - start_search - preproc_time
    return len(res), search_time, preproc_time

# T = "bbbbbbbbababbbbbbbb"
# P = "aba"
# res, search_time, preproc_time = search(P, T)
# for i in range(len(res)):
#     print(res[i], end=" ")
# print(f"Preprocessing time: {preproc_time:.8f}s")
# print(f"Search time: {search_time:.8f}s")

parser = argparse.ArgumentParser(
    description="Knuth Morris Pratt string matching algorithm"
)
parser.add_argument(
    "--input-file", "--sequence-file",
    dest="input_file",
    required=True,
    help="Path to input file (.txt or .fna)"
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

if args.pattern_file:
    with open(args.pattern_file, encoding="utf-8") as f:
        pattern = f.read().strip()
    if file_type == "fasta":
        pattern = pattern.replace("\n", "")
else:
    pattern = args.pattern.replace("\n", "\\n")

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
total_time = 0
for _ in range(args.repeat):
    match_count, matchtime, preproctime = search(pattern, text)
    total_time += matchtime+preproctime

print(f"Alphabet size:      {len(alphabet)}")
print(f"Preprocessing time: {preproctime:.8f} s")
print(f"Matching time:      {(matchtime):.8f} s")
print(f"Total time:         {(preproctime + matchtime):.8f} s")
print(f"Average time:       {(total_time/args.repeat):.8f} s")
print(f"Match count:        {match_count}")