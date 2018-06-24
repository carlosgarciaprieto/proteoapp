#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#"""Python 3.6.3 under macOS High Sierra 10.13.1"""


#"""My app import modules"""

from flask import Flask, render_template 
import monoisotopic
import prot
import Bio
from Bio import ExPASy
from Bio import SwissProt
from flask import request
app = Flask(__name__)


#"""Mass module computes monoisotopic mass of a protein sequence"""

@app.route('/')
def main():
    return render_template('index.html')

@app.route('/mass.html')
def mass():
    return render_template('mass.html')

@app.route('/mass_results', methods=['POST']) 
def mass_results():
    input_seq = request.form['input_seq'] 
    input_seq=input_seq.upper().strip().replace(' ', '').replace('\n', '').replace('\r', '').split(",")
    mass = request.form['mass']
    output = monoisotopic.masa(input_seq)
    return render_template('mass_results.html', **locals())


#"""Prot module implements Biopython capabilities regarding protein function, gene ontology, pfam domains and PDB structure visualization"""
    
@app.route('/uniprot.html')
def uniprot():
    return render_template('uniprot.html')   

@app.route('/uniprot_results', methods=['POST']) 
def uniprot_results():
    query_proteins = request.form['query_proteins'] 
    query_proteins = query_proteins.upper().strip().replace(' ', '').replace('\n', '').replace('\r', '').split(",")
    uniprot = request.form['uniprot']
    if uniprot == "function":
        output_list = prot.protfunction(query_proteins)
    elif uniprot == "ontology":
        output_list = prot.geneontology(query_proteins)
    elif uniprot == "pfam":
        output = prot.pfamdomains(query_proteins)
    elif uniprot == "pdb":
        output = prot.pdbstructure(query_proteins)
    return render_template('uniprot_results.html', **locals())
    

# Start the app (let's keep debug=True during debugging)
if __name__ == '__main__':
    app.run(debug=True)
