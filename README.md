# Gene_Edit_1
Custom gene editing code version 1

gene_editing.sh is intended to run directly on paired fastq sequencing files. 

    pear is used to join the R1 and R2 files

    a custom python script is run to filter reads for quality

    needle from EMBOSS is run to compare reads to a given reference sequence

    some command line editing is run on the needle output(SAM)

after the commmands of gene_editing.sh, a python script is used to summarize the SAM files

easyGeneEdit.py inputFile.sam outputfile.txt ref_sequence
