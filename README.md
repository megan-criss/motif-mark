# motif-mark
A python script to visualize motifs on sequences in comparison to exons.


## Introduction

The purpose of this script is to plot protein binding motifs on an image of the exon and flanking introns.

### Required Modules

This script requires the installation of pycairo for the data visualization portion.

'''

conda install -c conda-forge pycairo

'''

or 

'''

pip install pycairo

'''

This script also requires the installation  of argparse and re(regex):

'''

pip install argparse

'''

'''

pip install regex

'''



### How to Run

In order to run this code use:

'''

motif-mark.py -f (fasta file input) -m (text file containing motifs, each on a new line [example below])

'''

For this program, exons in the fasta files should be capitalized  in order for the program to properly demarkate them (example below).


#### Example Fasta file & motif file:

Fasta file formatting:

'''

\>Gene 1

agatgtgactctatcgctagcatcgatcgactACTGACGATCTACGATCactgatcatgcatgctagcatgc

'''

Motif file formatting:

'''

ygcy 
ATYG 
YYYYYY 
yyGCaa 

'''

### Output

This program outputs an svg file that is scale to the length of the genes in the fasta file. This program can take up to 7 different motifs.

SVG output will be named the same as the input fasta file


#### Example:

input:

'''

./motif-mark.py -f Figure_1.fasta -m Fig_1_motifs.txt

'''

output:

'''

Figure_1.svg

'''
