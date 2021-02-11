#!/usr/bin/env python

import argparse
import re
import cairo

#argparse expressions to take in user input
def get_args():
    parser = argparse.ArgumentParser(description = " A program to find motifs within a fasta file")
    parser.add_argument("-f", "--file", help ="SAM file", required = True)
    parser.add_argument("-m", "--motifs", help = "file containing motifs that you are searching for", required = True)
    return parser.parse_args()
args = get_args()
fasta_file = args.file
motif_file = args.motifs 

#this line determines output file name based on input file name
output_file = "".join(fasta_file.split('.')[:-1]) + '.svg'


#this defines a dictionary full of IUPAC degeneragte bases to be mathed to giiven motifs
# keys are IUPAC Bases, values are Regex expressions
iupac = {
        "A":"[Aa]",  #adenine
        "C":"[Cc]",  #cytosine
        "G":"[Gg]",  #Guanine
        "T":'[TtUu]',  #thymine/uracil
        "U":'[UuTt]',      #uracil/thymine
        "W":'[AaTtUu]',  #weak
        "S":'[CcGg]',  #strong
        "M":'[AaCc]',  #amino
        "K":'[GgTtUu]',  #keto
        "R":'[AaGg]',  #Purine
        "Y":'[CcTtUu]',  #pyrimidine
        "B":'[CcGgTtUu]',  #not A
        "D":'[AaGgTtUu]',  #not C
        "H":'[AaCcTtUu]',  #not G
        "V":'[AaCcGg]',  #not T
        "N":'[AaCcGgTtUu]',  #Any one base
        "Z":'[]'}   #zero

#this code defnes a list of all motifs given, to later be used for color matching, iupac_finder, etc.
with open(motif_file, 'r') as motif_file:
    motif_list = []
    for motif in motif_file:
        motif = motif.strip()
        motif_list.append(motif)


def iupac_finder(motif_list):
    '''
    A function to define the IUPAC nucleotide base pairs and return the regex expression for a given motif 
    '''
    conv_list = {}
    for motif in motif_list:
        conv_motif = ''
        motif = motif.upper()
        for mote in motif:
            if mote in iupac:
                conv_motif += iupac[mote]
            else:
                print("ERROR, mystery iupac detected: is that really what you wanted to do?")
        conv_list[conv_motif] = len(motif)
    return(conv_list)


#this code creates a dictionary with keys as the header of the fasta(lines that start with a > and contain the gene name)
#and the values of the dictionary as thee entire fasta sequence given, to be used for indexing and motif placement
with open(fasta_file, "r") as fasta:
    fasta_dict = {}
    gene_names = []
    for line in fasta:
        line = line.strip()
        if line.startswith('>'):
            header = line
            gene_names.append(header)
            read = ''
        else:
            read += line
        i = {header:read}
        fasta_dict.update(i)


def find_exon(line):
    '''
A function to find the exon coordinates in a fasta file
    '''
    caps = "([AUTCG]+)"
    line = line.strip()
    #print(line)
    if line.startswith(">"):
        pass
    else:
        indicies_exon = [(m.start(0), m.end(1)-1) for m in re.finditer(caps, line)]
        return(indicies_exon)

#defining all dictionaries that will be used and liast that will be used
exon_dict = {}
motif_dict = {}
color_dict_trans_motif = {}     #color dictionaryfor IUPAC translated motifs
color_dict = {}                 #color dictionary for original motif provided - colors will correspond to translated motifs
info_list = []

# this code block fills lists and calls the defined functions above to fill the lists/dictionaries that are specified above
#it also sets corresponding colors for the svg outut
for header,sequence in fasta_dict.items():
    exons = find_exon(sequence)
    e = {header:exons}
    exon_dict.update(e)
    motif_regex = iupac_finder(motif_list)
    color_list = [(0.2,0.4,1), (1,0.6,0), (1,0,0.1), (1,0.3,1), (1,1,0.2),(0.1,1,0),(0.6,0.2,0.7)]
                # Blue          orange      red         pink      yellow    green       purple
                # I encoded for more than required...just in case!
    LN = 0
    line = 0
    for motif, length in motif_regex.items():
        indicies_motif = [(m.start(0)) for m in re.finditer(rf'(?=({motif}))', sequence, re.IGNORECASE)]
        info_list.append((header, motif, length, indicies_motif))
        color_dict_trans_motif[motif] = (color_list[LN])
        LN += 1
    for motif in motif_list:
        color_dict[motif] = (color_list[line])
        line += 1


def draw_svg():
    '''
    This function is responsible for draing the svg output of the input fasta and motif files
    '''
    yax = 0
    width, height = 1000,1000
    surface = cairo.SVGSurface(output_file, width, height)
    context = cairo.Context(surface)
    for header, sequence in fasta_dict.items():
        sequence = sequence.strip()
        yax += 150
        xax = 65
        context.set_source_rgba(0, 0, 0)
        context.set_line_width(2)
        context.move_to(xax,yax)
        context.line_to(len(sequence)+xax,yax)
        context.stroke()
        x = exon_dict[header][0][0]
        y = exon_dict[header][0][1]
        context.set_source_rgb(0, 0, 0)
        context.set_line_width(30)
        context.move_to(x+xax,yax)
        context.line_to(y+xax,yax)
        context.stroke()
        for item in info_list:
            if item[0] == header:
                for m, length in motif_regex.items():
                    if m == item[1]:
                        r,b,g = color_dict_trans_motif[m]
                        for index in item[3]:
                            line_width = 15
                            context.set_source_rgb(r,b,g)
                            context.rectangle(index+xax, yax-7, length, line_width)
                            context.fill()
                            context.set_source_rgb(0,0,0)   
                            context.move_to(xax , yax -30 )
                            context.show_text(header)
    yleg = 700
    for motif in motif_list:
        x,y,z = color_dict[motif]
        #print(motif)
        context.set_source_rgb(0,0,0)
        context.move_to(65 , yleg+5)
        context.show_text(motif)
        context.set_source_rgb(x,y,z)
        context.move_to(50, yleg-2)
        context.line_to(60, yleg+2)
        context.stroke()
        yleg += 25
    context.set_source_rgb(0,0,0)
    context.move_to(65 , yleg+5)
    context.show_text("Exon")
    context.set_source_rgb(0,0,0)
    context.move_to(50, yleg-2)
    context.line_to(60, yleg+2)
    context.stroke()
    surface.finish()

#drawing the final output!!
draw_svg()
