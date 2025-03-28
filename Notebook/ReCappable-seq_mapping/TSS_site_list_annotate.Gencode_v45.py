import argparse
import pandas as pd
import scipy.stats
import numpy as np
from Bio import SeqIO
from Bio.Seq import reverse_complement
import regex
import warnings
import os
import time
from time import gmtime, strftime
warnings.filterwarnings("ignore")

bedtools = "/home/fox/Software/bedtools/2.29.1/bedtools"
blastn = "/home/fox/Software/blast/2.13.0/bin/blastn"

def mark_low_complexity(row):
    chr = str(row.name[0])
    pos = int(row.name[1])
    strand = row.name[2]
    
    reject = False

    if reference.get(chr) is None:
        return None
    if strand == "+":
        flanking = reference[chr][pos-3: pos+4].replace("A", "G")
        res = regex.finditer("(GGGGGGG){e<1}", flanking, overlapped=True)
        for i in res:
            reject = True
            break
        flanking = reference[chr][pos-1: pos+6].replace("A", "G")
        res = regex.finditer("(GGGGGGG){e<1}", flanking, overlapped=True)
        for i in res:
            reject = True
            break
        return reject
    else:
        flanking = reverse_complement(reference[chr][pos-3: pos+4]).replace("A", "G")
        res = regex.finditer("(GGGGGGG){e<1}", flanking, overlapped=True)            
        for i in res:
            reject = True
            break
        flanking = reference[chr][pos-7: pos].replace("A", "G")
        res = regex.finditer("(GGGGGGG){e<1}", flanking, overlapped=True)
        for i in res:
            reject = True
            break
        return reject
    
if __name__ == "__main__":
    description = """"""
    parser = argparse.ArgumentParser(prog="",fromfile_prefix_chars='@',description=description,formatter_class=argparse.RawTextHelpFormatter)
    #Require
    group_required = parser.add_argument_group("Required")
    group_required.add_argument("-i", dest="input", required=True, help="Input TSS BED files")
    group_required.add_argument("-o","--output",dest="output",required=False, default="TSS.out.bed",help="Output")
    group_required.add_argument("-f","--fasta",dest="fasta",required=False, default="/home/fox/Database/Genome/Human/Gencode/v45/GRCh38.primary_assembly.genome.nochr.fa")
    group_required.add_argument("-t","--tss",dest="tss_db",required=False, default="/home/fox/Database/Genome/Human/Gencode/v45/gencode.v45.primary_assembly.annotation.nochr.tss.bed",help="")
    group_required.add_argument("-e","--exons",dest="exons_db",required=False, default="/home/fox/Database/Genome/Human/Gencode/v45/gencode.v45.primary_assembly.annotation.nochr.tx.bed",help="")
    group_required.add_argument("-a","--annot",dest="annot",required=False, default="/home/fox/Database/Genome/Human/Gencode/v45/gencode.v45.primary_assembly.annotation.nochr.anno",help="")
    group_required.add_argument("-b","--blastdb",dest="blast_db",required=False, default="/home/fox/Database/BLAST/Human_snRNA_snoRNA/snRNA_snoRNA_blast",help="")
    group_required.add_argument("-r","--repeats",dest="repeats",required=False, default="/home/fox/Database/Repeats/hg38_repeats.nochr.sorted.bed",help="")
    options = parser.parse_args()

    print("[{t}] Analysis begins.".format(t=strftime("%Y-%m-%d %H:%M:%S", time.localtime())))
    
    reference = {}
    for seq in SeqIO.parse(options.fasta, "fasta"):
        reference[seq.id] = str(seq.seq).upper()

    print("[{t}] Reference loaded.".format(t=strftime("%Y-%m-%d %H:%M:%S", time.localtime())))
    
    data = {}
    all_sites = {}

    print("[{t}] Searching for snRNA/snoRNA like sequences...".format(t=strftime("%Y-%m-%d %H:%M:%S", time.localtime())))
    
    df = pd.read_csv(options.input, sep="\t", index_col=None, names=["chr", "pos_0", "pos_1", "x", "xx", "strand"], low_memory=False)

    with open(options.output + ".site", "w") as output:
        for idx, row in df.iterrows():
            output.write("{}\t{}\t{}\n".format(row["chr"], row["pos_1"], row["strand"]))
    
    print("[{t}] Getting fasta...".format(t=strftime("%Y-%m-%d %H:%M:%S", time.localtime())))
    cmd = "/home/fox/Software/bin/python /home/fox/Scripts/genome_flanking_v3.py -U 0 -D 50 -s Gencodev45 -i {sites_in} > {fasta_out}".format(
        sites_in=options.output + ".site",
        fasta_out=options.output + ".fa")
    print(cmd)
    # os.system(cmd)
    
    df = df.drop(["x", "xx"], axis=1)
    df["chr"] = df["chr"].astype(str)
    df = df.set_index(["chr", "pos_1", "strand"])
   
    # os.system("{blastn} -db {blast_db} -out {blast_out} -outfmt 6 -query {fasta_in} -num_threads 10 -qcov_hsp_perc 50 -perc_identity 50 -word_size 10".format(blastn=blastn, blast_db=options.blast_db, blast_out=options.output+".blast_out.txt", fasta_in=options.output + ".fa"))

    snRNA_snoRNA_keys = []
    with open(options.output+".blast_out.txt", "r") as input:
        for line in input.readlines():
            line = line.strip().split("\t")
            if float(line[-1]) >= 40:
                c, p, s = line[0].split("@")
                snRNA_snoRNA_keys.append((c, int(p), s))
    
    df["snRNA_snoRNA_like"] = False
    df.loc[snRNA_snoRNA_keys, "snRNA_snoRNA_like"] = True

    # Annotate to closest TSS
    # Extract sites to BED file first

    print("[{t}] Proceed to gene annotation...".format(t=strftime("%Y-%m-%d %H:%M:%S", time.localtime())))
    
    with open(options.output + ".bed", "w") as bed_out:
        for idx, row in df.iterrows():
            bed_out.write("{}\t{}\t{}\t{}\t{}\t{}\n".format(idx[0], idx[1] - 1, idx[1], ".", ".", idx[2]))

    # Bedtools closest
    cmd = "{bedtools} closest -s -D b -a {bed_out} -b {database} > {closest_tss}".format(bedtools=bedtools, bed_out = options.input, closest_tss= options.output + ".closest_tss.bed", database=options.tss_db)
    print(cmd)
    # os.system(cmd)
    

    cmd = "{bedtools} closest -S -D b -a {bed_out} -b {database} > {closest_tss}".format(bedtools=bedtools, bed_out = options.input, closest_tss= options.output + ".closest_tss.upstream.bed", database=options.tss_db)
    print(cmd)
    # os.system(cmd)
    

    cmd = "{bedtools} closest -s -D b -a {bed_out} -b {database} > {closest_exons}".format(bedtools=bedtools, bed_out =options.input, closest_exons= options.output + ".closest_exons.bed", database=options.exons_db)
    print(cmd)
    # os.system(cmd)
    

    # Merge annotations
    print("[{t}] Collasping annotations...".format(t=strftime("%Y-%m-%d %H:%M:%S", time.localtime())))
    cmd = "python /home/fox/Scripts/ReCappable-seq/collaspe_bed_annotations_v3.py -i {closest} -a {annot} -o {collasped}".format(closest= options.output + ".closest_tss.bed", annot=options.annot, collasped = options.output + ".collasped.1.bed")
    print(cmd)
    # os.system(cmd)

    cmd = "python /home/fox/Scripts/ReCappable-seq/collaspe_bed_annotations_fix_other_exons_v2.py -a {annot} -i {collasped} -I {closest_exons} -o {collasped_2} ".format(collasped = options.output + ".collasped.1.bed", closest_exons=options.output + ".closest_exons.bed", annot=options.annot, collasped_2= options.output + ".collasped.2.bed")
    print(cmd)
    os.system(cmd)

    cmd = "python /home/fox/Scripts/CROWN-seq/collaspe_bed_annotations_fix_upstream_TSS.py  --thresh-min -2000 --thresh-max 2000 -a {annot} -i {collasped} -I {closest_exons} -o {collasped_final} ".format(collasped = options.output + ".collasped.2.bed", closest_exons=options.output + ".closest_tss.upstream.bed", annot=options.annot, collasped_final= options.output + ".collasped_final.bed")
    print(cmd)
    os.system(cmd)

    # Repeat mask
    cmd = "{bedtools} intersect -wa -wb -a {bed_in} -b {bed_repeats} > {bed_out}".format(bedtools = bedtools, bed_in = options.input, bed_out = options.output + ".repeats.bed", bed_repeats = options.repeats)
    print(cmd)
    os.system(cmd) 


    data_annot_gene = {}
    data_annot_biotype = {}
    data_annot_method = {}
    data_annot_dist = {}
    data_annot_repeat = {}

    with open(options.output + ".collasped_final.bed", "r") as annotated_bed:
        for line in annotated_bed.readlines():
            line = line.strip().split("\t")
            key = (line[0], int(line[2]), line[5])
            annot = line[3]
            gene, gene_biotype, trans_biotype, method, dist = line[3].split("|")

            data_annot_gene[key] = gene
            data_annot_biotype[key] = gene_biotype
            data_annot_method[key] = method
            data_annot_dist[key] = dist

    with open(options.output + ".repeats.bed", "r") as annotated_bed:
        for line in annotated_bed.readlines():
            line = line.strip().split("\t")
            key = (line[0], int(line[2]), line[5])
            annot = line[-2]
            data_annot_repeat[key] = annot

    
    # print(df)
    df["Gene"] = df.apply(lambda x: data_annot_gene.get((x.name[0], x.name[1], x.name[2])), axis=1)
    df["Gene_biotype"] = df.apply(lambda x: data_annot_biotype.get((x.name[0], x.name[1], x.name[2])), axis=1)
    df["Annotate_to_closest"] = df.apply(lambda x: data_annot_method.get((x.name[0], x.name[1], x.name[2])), axis=1)
    df["Relative_distance"] = df.apply(lambda x: data_annot_dist.get((x.name[0], x.name[1], x.name[2])), axis=1)
    df["Repeat"] = df.apply(lambda x: data_annot_repeat.get((x.name[0], x.name[1], x.name[2])), axis=1)

    df.to_csv(options.output)
    # df_temp.to_csv(options.output + ".test.csv")
    print("[{t}] DONE!".format(t=strftime("%Y-%m-%d %H:%M:%S", time.localtime())))



