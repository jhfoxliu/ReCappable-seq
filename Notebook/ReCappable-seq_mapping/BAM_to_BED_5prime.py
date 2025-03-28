import pysam
import sys
from collections import defaultdict

# Single end
data = defaultdict(int)
N = 0
with pysam.AlignmentFile(sys.argv[1], "rb") as SAM:
    for read in SAM:
        if read.is_secondary == False and read.is_unmapped == False:
            if read.is_read1 and read.is_reverse == False:
                align_pairs = read.get_aligned_pairs(matches_only=False)
                i = 0
                while True:
                    useful_base = align_pairs[i]
                    if useful_base[1] is not None:
                        key = (read.reference_name, useful_base[1], "+") # 0-based 
                        data[key] += 1
                        N += 1
                        # sys.stderr.write("{}\r".format(N))
                        break
                    i += 1
            elif read.is_read1 and read.is_reverse == True:
                align_pairs = read.get_aligned_pairs(matches_only=False)
                i = -1
                while True:
                    useful_base = align_pairs[i]
                    if useful_base[1] is not None:
                        key = (read.reference_name, useful_base[1], "-") # 0-based 
                        data[key] += 1
                        N += 1
                        break
                        # sys.stderr.write("{}\r".format(N))
                    i -= 1
for key, values in data.items():
    print("{}\t{}\t{}\t{}\t{}\t{}".format(key[0], key[1], key[1] + 1, ".", values, key[2]))
