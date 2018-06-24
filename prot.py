#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#"""Python 3.6.3 under macOS High Sierra 10.13.1"""

"""prot module 
   contains functions for implementing some Biopython capabilities regarding protein function, gene ontology, Pfam domains and PDB structure visualization"""

#"""We import Biopython module"""

import Bio
from Bio import ExPASy
from Bio import SwissProt


#"""We define the function to show the protein function information"""

def protfunction(query_proteins):
    """Shows the proteins function given their names or ids
    str -> list"""
    function_list = []
    for prot in query_proteins:
        with ExPASy.get_sprot_raw(prot) as handle:
            record = SwissProt.read(handle)
            function_list.append((prot,record.comments[0][10:]))
    return function_list


#"""We define the function to retrieve the biological processes related to proteins"""

def geneontology(query_proteins):
    """Retrieves gene ontology biological processes given protein names or ids
    str -> set"""
    gene_ontology = []
    for prot in query_proteins:
        with ExPASy.get_sprot_raw(prot) as handle:
            record = SwissProt.read(handle)
            for ref in record.cross_references:
                if ref[0]=="GO" and ref[2].startswith("P"):
                    gene_ontology.append((prot,ref[2].split(":")[1]))
    return gene_ontology


#"""We define the function to retrieve the Pfam domains of the proteins"""

def pfamdomains(query_proteins):
    """Retrieves pfam domains given protein names or ids
    str -> list"""
    global pfam_list
    pfam_list = []
    template_url="https://pfam.xfam.org/family/{}"
    for prot in query_proteins:
        with ExPASy.get_sprot_raw(prot) as handle:
            record = SwissProt.read(handle)
            for ref in record.cross_references:
                if ref[0]=="Pfam":
                    pfam_list.append((prot,ref[2],template_url.format(ref[2])))
    return pfam_list


#"""We define the function to retrieve the PDB structure information of the proteins"""

def pdbstructure(query_proteins):
    """Retrieves PDB structure information given protein names or ids
    str -> list"""
    global pdb_list
    pdb_list = []
    template_url="https://www.rcsb.org/structure/{}"
    template_url_images="https://cdn.rcsb.org/images/rutgers/{2}/{1}/{0}.pdb1-500.jpg"
    for prot in query_proteins:
        with ExPASy.get_sprot_raw(prot) as handle:
            record = SwissProt.read(handle)
            for ref in record.cross_references:
                if ref[0]=="PDB":
                    pdb_list.append((prot,ref[1:4],template_url.format(ref[1]),template_url_images.format(ref[1].lower(),ref[1].lower(),(ref[1].lower())[1:-1])))
    return pdb_list

 
