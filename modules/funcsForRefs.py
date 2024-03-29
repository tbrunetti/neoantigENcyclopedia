import sys
import os
sys.path.append(os.getcwd())
from typing import *
from Codon import *
from NucleicAcids import *

def generate_codon_reference() -> Dict[str, Codon]:
    # generate amino acid library
    codon_library = {} 
    
    # reference information for all 64 possible codons       
    codon_library['UUU'] = Codon(fullname = 'Phenylalanine', shortname = 'Phe', symbol ='F', hydropathy = 'hydrophobic', chemical_class = 'aromatic', charge = 'uncharged', hydrogen_donor_status = 'none', polarity = 'nonpolar')
    codon_library['UUC'] = Codon(fullname = 'Phenylalanine', shortname = 'Phe', symbol ='F', hydropathy = 'hydrophobic', chemical_class = 'aromatic', charge = 'uncharged', hydrogen_donor_status = 'none', polarity = 'nonpolar')
    codon_library['UUA'] = Codon(fullname = 'Leucine', shortname = 'Leu', symbol = 'L', hydropathy = 'hydrophobic', chemical_class = 'aliphatic', charge = 'uncharged', hydrogen_donor_status = 'none', polarity = 'polar')
    codon_library['UUG'] = Codon(fullname = 'Leucine', shortname = 'Leu', symbol = 'L', hydropathy = 'hydrophobic', chemical_class = 'aliphatic', charge = 'uncharged', hydrogen_donor_status = 'none', polarity = 'polar')
    codon_library['UCU'] = Codon(fullname = 'Serine',shortname = 'Ser', symbol = 'S', hydropathy = 'neutral', chemical_class= 'hydroxyl', charge = 'uncharged', hydrogen_donor_status = 'donor and acceptor', polarity = 'polar')
    codon_library['UCC'] = Codon(fullname = 'Serine',shortname = 'Ser', symbol = 'S', hydropathy = 'neutral', chemical_class= 'hydroxyl', charge = 'uncharged', hydrogen_donor_status = 'donor and acceptor', polarity = 'polar')
    codon_library['UCA'] = Codon(fullname = 'Serine',shortname = 'Ser', symbol = 'S', hydropathy = 'neutral', chemical_class= 'hydroxyl', charge = 'uncharged', hydrogen_donor_status = 'donor and acceptor', polarity = 'polar')
    codon_library['UCG'] = Codon(fullname = 'Serine',shortname = 'Ser', symbol = 'S', hydropathy = 'neutral', chemical_class= 'hydroxyl', charge = 'uncharged', hydrogen_donor_status = 'donor and acceptor', polarity = 'polar')
    codon_library['UAU'] = Codon(fullname = 'Tyrosine', shortname = 'Tyr', symbol = 'Y', hydropathy = 'neutral', chemical_class = 'aromatic', charge = 'uncharged', hydrogen_donor_status = 'donor and acceptor',polarity = 'polar')
    codon_library['UAC'] = Codon(fullname = 'Tyrosine', shortname = 'Tyr', symbol = 'Y', hydropathy = 'neutral', chemical_class = 'aromatic', charge = 'uncharged', hydrogen_donor_status = 'donor and acceptor',polarity = 'polar')
    codon_library['UAA'] = Codon(fullname = 'STOP', shortname = 'STOP', symbol = '*', hydropathy = 'n/a', chemical_class = 'n/a', charge = 'n/a', hydrogen_donor_status = 'n/a', polarity = 'n/a')
    codon_library['UAG'] = Codon(fullname = 'STOP', shortname = 'STOP', symbol = '*', hydropathy = 'n/a', chemical_class = 'n/a', charge = 'n/a', hydrogen_donor_status = 'n/a', polarity = 'n/a')
    codon_library['UGU'] = Codon(fullname = 'Cysteine', shortname = 'Cys', symbol = 'C', hydropathy = 'hydrophobic', chemical_class = 'sulfur', charge = 'uncharged', hydrogen_donor_status = 'none', polarity = 'nonpolar')
    codon_library['UGC'] = Codon(fullname = 'Cysteine', shortname = 'Cys', symbol = 'C', hydropathy = 'hydrophobic', chemical_class = 'sulfur', charge = 'uncharged', hydrogen_donor_status = 'none', polarity = 'nonpolar')
    codon_library['UGA'] = Codon(fullname = 'STOP', shortname = 'STOP', symbol = '*', hydropathy = 'n/a', chemical_class = 'n/a', charge = 'n/a', hydrogen_donor_status = 'n/a', polarity = 'n/a')
    codon_library['UGG'] = Codon(fullname = 'Tryptophan', shortname = 'Trp', symbol = 'W', hydropathy = 'hydrophobic', chemical_class = 'aromatic', charge = 'uncharged', hydrogen_donor_status = 'donor', polarity = 'nonpolar')
    codon_library['CUU'] = Codon(fullname = 'Leucine', shortname = 'Leu', symbol = 'L', hydropathy = 'hydrophobic', chemical_class = 'aliphatic', charge = 'uncharged', hydrogen_donor_status = 'none', polarity = 'nonpolar')
    codon_library['CUC'] = Codon(fullname = 'Leucine', shortname = 'Leu', symbol = 'L', hydropathy = 'hydrophobic', chemical_class = 'aliphatic', charge = 'uncharged', hydrogen_donor_status = 'none', polarity = 'nonpolar')
    codon_library['CUA'] = Codon(fullname = 'Leucine', shortname = 'Leu', symbol = 'L', hydropathy = 'hydrophobic', chemical_class = 'aliphatic', charge = 'uncharged', hydrogen_donor_status = 'none', polarity = 'nonpolar')
    codon_library['CUG'] = Codon(fullname = 'Leucine', shortname = 'Leu', symbol = 'L', hydropathy = 'hydrophobic', chemical_class = 'aliphatic', charge = 'uncharged', hydrogen_donor_status = 'none', polarity = 'nonpolar')
    codon_library['CCU'] = Codon(fullname = 'Proline', shortname = 'Pro', symbol = 'P', hydropathy = 'neutral', chemical_class = 'aliphatic', charge = 'uncharged', hydrogen_donor_status = 'none', polarity = 'nonpolar')
    codon_library['CCC'] = Codon(fullname = 'Proline', shortname = 'Pro', symbol = 'P', hydropathy = 'neutral', chemical_class = 'aliphatic', charge = 'uncharged', hydrogen_donor_status = 'none', polarity = 'nonpolar')
    codon_library['CCA'] = Codon(fullname = 'Proline', shortname = 'Pro', symbol = 'P', hydropathy = 'neutral', chemical_class = 'aliphatic', charge = 'uncharged', hydrogen_donor_status = 'none', polarity = 'nonpolar')
    codon_library['CCG'] = Codon(fullname = 'Proline', shortname = 'Pro', symbol = 'P', hydropathy = 'neutral', chemical_class = 'aliphatic', charge = 'uncharged', hydrogen_donor_status = 'none', polarity = 'nonpolar')
    codon_library['CAU'] = Codon(fullname = 'Histidine', shortname = 'His', symbol = 'H', hydropathy = 'neutral', chemical_class = 'basic', charge = 'positive', hydrogen_donor_status = 'donor and acceptor', polarity = 'polar')
    codon_library['CAC'] = Codon(fullname = 'Histidine', shortname = 'His', symbol = 'H', hydropathy = 'neutral', chemical_class = 'basic', charge = 'positive', hydrogen_donor_status = 'donor and acceptor', polarity = 'polar')
    codon_library['CAA'] = Codon(fullname = 'Glutamine', shortname = 'Gln', symbol = 'Q', hydropathy = 'hydrophilic', chemical_class = 'amide', charge = 'uncharged', hydrogen_donor_status = 'donor and acceptor', polarity = 'polar')
    codon_library['CAG'] = Codon(fullname = 'Glutamine', shortname = 'Gln', symbol = 'Q', hydropathy = 'hydrophilic', chemical_class = 'amide', charge = 'uncharged', hydrogen_donor_status = 'donor and acceptor', polarity = 'polar')
    codon_library['CGU'] = Codon(fullname = 'Arginine', shortname = 'Arg', symbol = 'R', hydropathy = 'hydrophilic', chemical_class = 'basic', charge = 'positive', hydrogen_donor_status = 'donor', polarity = 'polar')
    codon_library['CGC'] = Codon(fullname = 'Arginine', shortname = 'Arg', symbol = 'R', hydropathy = 'hydrophilic', chemical_class = 'basic', charge = 'positive', hydrogen_donor_status = 'donor', polarity = 'polar')
    codon_library['CGA'] = Codon(fullname = 'Arginine', shortname = 'Arg', symbol = 'R', hydropathy = 'hydrophilic', chemical_class = 'basic', charge = 'positive', hydrogen_donor_status = 'donor', polarity = 'polar')
    codon_library['CGG'] = Codon(fullname = 'Arginine', shortname = 'Arg', symbol = 'R', hydropathy = 'hydrophilic', chemical_class = 'basic', charge = 'positive', hydrogen_donor_status = 'donor', polarity = 'polar')
    codon_library['AUU'] = Codon(fullname = 'Isoleucine', shortname = 'Ile', symbol = 'I', hydropathy = 'hydrophobic', chemical_class = 'aliphatic', charge = 'uncharged', hydrogen_donor_status = 'none', polarity = 'nonpolar')
    codon_library['AUC'] = Codon(fullname = 'Isoleucine', shortname = 'Ile', symbol = 'I', hydropathy = 'hydrophobic', chemical_class = 'aliphatic', charge = 'uncharged', hydrogen_donor_status = 'none', polarity = 'nonpolar')
    codon_library['AUA'] = Codon(fullname = 'Isoleucine', shortname = 'Ile', symbol = 'I', hydropathy = 'hydrophobic', chemical_class = 'aliphatic', charge = 'uncharged', hydrogen_donor_status = 'none', polarity = 'nonpolar')
    codon_library['AUG'] = Codon(fullname = 'Methionine', shortname = 'Met', symbol = 'M', hydropathy = 'hydrophobic', chemical_class = 'sulfur', charge = 'uncharged', hydrogen_donor_status = 'none', polarity = 'nonpolar')
    codon_library['ACU'] = Codon(fullname = 'Threonine', shortname = 'Thr', symbol = 'T', hydropathy = 'neutral', chemical_class = 'hydroxyl', charge = 'uncharged', hydrogen_donor_status = 'donor and acceptor', polarity = 'polar')
    codon_library['ACC'] = Codon(fullname = 'Threonine', shortname = 'Thr', symbol = 'T', hydropathy = 'neutral', chemical_class = 'hydroxyl', charge = 'uncharged', hydrogen_donor_status = 'donor and acceptor', polarity = 'polar')
    codon_library['ACA'] = Codon(fullname = 'Threonine', shortname = 'Thr', symbol = 'T', hydropathy = 'neutral', chemical_class = 'hydroxyl', charge = 'uncharged', hydrogen_donor_status = 'donor and acceptor', polarity = 'polar')
    codon_library['ACG'] = Codon(fullname = 'Threonine', shortname = 'Thr', symbol = 'T', hydropathy = 'neutral', chemical_class = 'hydroxyl', charge = 'uncharged', hydrogen_donor_status = 'donor and acceptor', polarity = 'polar')
    codon_library['AAU'] = Codon(fullname = 'Asparagine', shortname = 'Asn', symbol = 'N', hydropathy = 'hydrophilic', chemical_class = 'amide', charge = 'uncharged', hydrogen_donor_status = 'donor and acceptor', polarity = 'polar')
    codon_library['AAC'] = Codon(fullname = 'Asparagine', shortname = 'Asn', symbol = 'N', hydropathy = 'hydrophilic', chemical_class = 'amide', charge = 'uncharged', hydrogen_donor_status = 'donor and acceptor', polarity = 'polar')
    codon_library['AAA'] = Codon(fullname = 'Lysine', shortname = 'Lys', symbol = 'K', hydropathy = 'hydrophilic', chemical_class = 'basic', charge = 'positive', hydrogen_donor_status = 'donor', polarity = 'polar')
    codon_library['AAG'] = Codon(fullname = 'Lysine', shortname = 'Lys', symbol = 'K', hydropathy = 'hydrophilic', chemical_class = 'basic', charge = 'positive', hydrogen_donor_status = 'donor', polarity = 'polar')
    codon_library['AGU'] = Codon(fullname = 'Serine',shortname = 'Ser', symbol = 'S', hydropathy = 'neutral', chemical_class= 'hydroxyl', charge = 'uncharged', hydrogen_donor_status = 'donor and acceptor', polarity = 'polar')
    codon_library['AGC'] = Codon(fullname = 'Serine',shortname = 'Ser', symbol = 'S', hydropathy = 'neutral', chemical_class= 'hydroxyl', charge = 'uncharged', hydrogen_donor_status = 'donor and acceptor', polarity = 'polar')
    codon_library['AGA'] = Codon(fullname = 'Arginine', shortname = 'Arg', symbol = 'R', hydropathy = 'hydrophilic', chemical_class = 'basic', charge = 'positive', hydrogen_donor_status = 'donor', polarity = 'polar')
    codon_library['AGG'] = Codon(fullname = 'Arginine', shortname = 'Arg', symbol = 'R', hydropathy = 'hydrophilic', chemical_class = 'basic', charge = 'positive', hydrogen_donor_status = 'donor', polarity = 'polar')
    codon_library['GUU'] = Codon(fullname = 'Valine', shortname = 'V', symbol = 'V', hydropathy = 'hydrophobic', chemical_class = 'aliphatic', charge = 'uncharged', hydrogen_donor_status = 'none', polarity = 'nonpolar')
    codon_library['GUC'] = Codon(fullname = 'Valine', shortname = 'V', symbol = 'V', hydropathy = 'hydrophobic', chemical_class = 'aliphatic', charge = 'uncharged', hydrogen_donor_status = 'none', polarity = 'nonpolar')
    codon_library['GUA'] = Codon(fullname = 'Valine', shortname = 'V', symbol = 'V', hydropathy = 'hydrophobic', chemical_class = 'aliphatic', charge = 'uncharged', hydrogen_donor_status = 'none', polarity = 'nonpolar')
    codon_library['GUG'] = Codon(fullname = 'Valine', shortname = 'V', symbol = 'V', hydropathy = 'hydrophobic', chemical_class = 'aliphatic', charge = 'uncharged', hydrogen_donor_status = 'none', polarity = 'nonpolar')
    codon_library['GCU'] = Codon(fullname = 'Alanine', shortname = 'Ala', symbol = 'A', hydropathy = 'hydrophobic', chemical_class = 'aliphatic', charge = 'uncharged', hydrogen_donor_status = 'none', polarity = 'nonpolar')
    codon_library['GCC'] = Codon(fullname = 'Alanine', shortname = 'Ala', symbol = 'A', hydropathy = 'hydrophobic', chemical_class = 'aliphatic', charge = 'uncharged', hydrogen_donor_status = 'none', polarity = 'nonpolar')
    codon_library['GCA'] = Codon(fullname = 'Alanine', shortname = 'Ala', symbol = 'A', hydropathy = 'hydrophobic', chemical_class = 'aliphatic', charge = 'uncharged', hydrogen_donor_status = 'none', polarity = 'nonpolar')
    codon_library['GCG'] = Codon(fullname = 'Alanine', shortname = 'Ala', symbol = 'A', hydropathy = 'hydrophobic', chemical_class = 'aliphatic', charge = 'uncharged', hydrogen_donor_status = 'none', polarity = 'nonpolar')
    codon_library['GAU'] = Codon(fullname = 'Aspartic acid', shortname = 'Asp', symbol = 'D', hydropathy = 'hydrophilic', chemical_class = 'acidic', charge = 'negative', hydrogen_donor_status = 'acceptor', polarity = 'polar')
    codon_library['GAC'] = Codon(fullname = 'Aspartic acid', shortname = 'Asp', symbol = 'D', hydropathy = 'hydrophilic', chemical_class = 'acidic', charge = 'negative', hydrogen_donor_status = 'acceptor', polarity = 'polar')
    codon_library['GAA'] = Codon(fullname = 'Glutamic acid', shortname = 'Glu', symbol = 'E', hydropathy = 'hydrophilic', chemical_class = 'acidic', charge = 'negative', hydrogen_donor_status= 'acceptor', polarity = 'polar')
    codon_library['GAG'] = Codon(fullname = 'Glutamic acid', shortname = 'Glu', symbol = 'E', hydropathy = 'hydrophilic', chemical_class = 'acidic', charge = 'negative', hydrogen_donor_status= 'acceptor', polarity = 'polar')
    codon_library['GGU'] = Codon(fullname = 'Glycine', shortname = 'Gly', symbol = 'G', hydropathy = 'neutral', chemical_class = 'aliphatic', charge = 'uncharged', hydrogen_donor_status = 'none', polarity = 'nonpolar')
    codon_library['GGC'] = Codon(fullname = 'Glycine', shortname = 'Gly', symbol = 'G', hydropathy = 'neutral', chemical_class = 'aliphatic', charge = 'uncharged', hydrogen_donor_status = 'none', polarity = 'nonpolar')
    codon_library['GGA'] = Codon(fullname = 'Glycine', shortname = 'Gly', symbol = 'G', hydropathy = 'neutral', chemical_class = 'aliphatic', charge = 'uncharged', hydrogen_donor_status = 'none', polarity = 'nonpolar')
    codon_library['GGG'] = Codon(fullname = 'Glycine', shortname = 'Gly', symbol = 'G', hydropathy = 'neutral', chemical_class = 'aliphatic', charge = 'uncharged', hydrogen_donor_status = 'none', polarity = 'nonpolar')
    
    return codon_library

def generate_fasta_sequences(fasta_file : str) -> Dict[str, Dna]:
    dna_fasta = {}
    
    with open(fasta_file, 'r') as reference:
        chromosome = ''
        for line in reference:
            if line.startswith('>'):
                chromosome = line.split()[0][1:]
                dna_fasta[chromosome] = ''
            else:
                dna_fasta[chromosome] = dna_fasta[chromosome] + line.strip().upper()
    
    return dna_fasta
