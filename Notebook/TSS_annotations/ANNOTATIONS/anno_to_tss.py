import sys

with open(sys.argv[1], "r") as input:
    for line in input.readlines():
        line = line.strip().split("\t")
        enst = line[1]
        chr = line[2]
        strand = line[3]
        txStart = int(line[4])
        txEnd = int(line[5])
        gene = line[-3]
        if chr == "chrM":
            chr = "MT"
        elif chr.startswith("chr"):
            chr = chr.replace("chr", "")
        if strand == "+":
            print("{}\t{}\t{}\t{}\t{}\t{}".format(chr, txStart, txStart+1,  gene, enst,strand))
        elif strand == "-":
            print("{}\t{}\t{}\t{}\t{}\t{}".format(chr, txEnd-1, txEnd,  gene, enst,strand))
    
