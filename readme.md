### About
This script takes as input a text file of DNA sequences separated by new lines. The script will output up to 4
pieces of information for each input sequence:
1) The initial sequence, the GC content of the initial sequence, and the translation of the original sequence
2) A shuffled version of the sequence where each codon is replaced with the most frequently used codon (as reported in https://www.kazusa.or.jp/codon/cgi-bin/showcodon.cgi?species=7227), the new shuffled sequences GC content, and a translation of the new sequence
3) A randomized version of the sequence, its GC content and its translation
4) A shuffled version where each codon is replaced with the codon for the same amino acid that most closely matches the original codons GC content.

### Input arguments
The script takes 4 input arguments:
1) The input file
2) The directory to save the output file (-o, --out_dir, default = /out)
3) The minimum acceptable GC content difference between the input and output sequence as a percent (-g, --gc_difference, default = 5). If the codon optimized version meets this cut-off outputs 3 and 4 will not be produced
4) The number of iterations for the randomly generated sequence (-i, --iterations, default = 1000)

### To do:
Add test suite
