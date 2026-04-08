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

def rabin_karp(file, pattern, d, q):
    text = loadFile(file)
    print("Text: ", file)
    print("Number of words:", len(text.split()))
    print("Number of unique words:", len(set(text.split())))
    n = len(text)
    m = len(pattern)
    h = pow(d, m-1, q)
    p = 0
    ts = 0
    validShifts = 0
    for i in range(m):
        p = (d*p+ ord(pattern[i])) % q
        ts = (d*ts+ ord(text[i])) % q
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
    print(f"Valid shifts: {validShifts}")



for item in [
    ["text1.txt", "AW", 256, 101],
    ["dna1000.txt", "gaa", 4, 101],
    ["newsarticle.txt", "help", 256, 9973]
]:
    start = time.perf_counter()
    rabin_karp(item[0], item[1], item[2], item[3])
    print(f"Time taken: {(time.perf_counter()-start):.9f}\n")