# ReCappable-seq

The ReCappable-seq analysis pipeline includes (1) the mapping of ReCappable-seq reads; (2) annotate the TSSs to the closest known TSS.

(1) Example ReCappable-seq mapping pipeline can be found under `Notebook/RecCappable-seq_mapping/ReCappable-seq_mapping.ipynb`. 

(2) Example TSS annotation pipeline can be found under `Notebook/TSS_annotations/TSS_annotations.ipynb`.

The step (2) will end up with a CSV file containing the annotations of the TSS found by the mapping pipeline. You can merge the coverage from step (1) by the table.

(please unzip `hg38_repeats.nochr.sorted.bed.gz`, `chr1.fa.gz`, and `chr1.gtf.gz` first in the `ANNOTATIONS` folder first)



