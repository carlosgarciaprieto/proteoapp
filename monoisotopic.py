#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#"""Python 3.6.3 under macOS High Sierra 10.13.1"""

"""monoisotopic module 
contains masa function to compute the monoisotopic mass of a given protein sequence"""

#"""We define masa function to compute the monoisotopic mass of a protein sequence"""

def masa(input_seq):
    """Computes the monoisotopic mass (Da) of proteins given their sequences
    str -> str"""
    output_list=[]
    monoisotopic_mass = dict([["A", 71.03711],
                              ["C",   103.00919],
                              ["D",   115.02694],
                              ["E",   129.04259],
                              ["F",   147.06841],
                              ["G",    57.02146],
                              ["H",   137.05891],
                              ["I",   113.08406],
                              ["K",   128.09496],
                              ["L",   113.08406],
                              ["M",   131.04049],
                              ["N",   114.04293],
                              ["P",    97.05276],
                              ["Q",   128.05858],
                              ["R",   156.10111],
                              ["S",    87.03203],
                              ["T",   101.04768],
                              ["V",    99.06841],
                              ["W",   186.07931],
                              ["Y",   163.06333]]) 
    for protein in input_seq:
        weight = sum(monoisotopic_mass[aa] for aa in protein)
        output_text = "The monoisotopic mass of the protein sequence {} is {} Da".format(protein,weight)
        output_list.append(output_text) 
    return output_list
