import time
import os

def loadFile(filename):
    text = ""
    # Get this directory
    base_dir = os.path.dirname(os.path.abspath(__file__))
    # Navigate up one level, then into Texts
    texts_dir = os.path.join(base_dir, "..", "Texts")
    with open(os.path.join(texts_dir, filename), encoding="utf-8") as f:
        for line in f:
            text += line.strip()
    return text

def naive(file, pattern):
    text = loadFile(file)
    print("Text: ", file)
    print("Number of words:", len(text.split()))
    print("Number of unique words:", len(set(text.split())))
    n = len(text)
    m = len(pattern)
    validshifts = 0
    for s in range(0, n-m+1):
        patternMatch = True
        for i in range(m):
            if(text[s+i] != pattern[i]):
                patternMatch = False
                break
        if(patternMatch):
            validshifts += 1
            #print(f"Valid shift {s}")
    print(f"Valid shifts: {validshifts}")

for item in [
    ["text1.txt", "AW"],
    ["dna1000.txt", "gaa"],
    ["newsarticle.txt", "help"]
]:
    start = time.perf_counter()
    naive(item[0], item[1])
    print(f"Time taken: {(time.perf_counter()-start):.9f}\n")