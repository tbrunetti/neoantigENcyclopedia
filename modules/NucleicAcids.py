from typing import *
import sys
import os
sys.path.append(os.getcwd())
from Codon import *


class Rna:
    
    rna_pairs = {'A':'U', 'U':'A', 'G':'C', 'C':'G'}
    
    RNA_NUCLEOTIDES = ('A', 'U', 'G', 'C')
    
    def __init__(self, sequence : str, rna_base_options = RNA_NUCLEOTIDES) -> None:
        check_nt = set([sequence[pos].upper() for pos in range(0, len(sequence))])
        if check_nt <= set(rna_base_options):
            print('Confirmed sequence contains exclusive RNA basepairs')
        else:
            raise ValueError('%s is not a valid nucleotide for a RNA object.' % check_nt)
        
        self.sequence = sequence.upper()
        
        def translate(self, codon_library : Dict[str, Codon]) -> str:
            protein_seq = []
            for codon in range(0, len(self.sequence), 3):
                print(codon)
                protein_seq.append(codon_library[self.sequence[codon:codon+3]].translate_symbol())
            
            return ''.join(protein_seq)
        
            
class Dna:
    dna_pairs = {'A':'T', 'T':'A', 'G':'C', 'C':'G'}
    
    DNA_NUCLEOTIDES = ('A', 'T', 'G', 'C')
    
    def __init__(self, sequence : str, dna_base_options = DNA_NUCLEOTIDES) -> None:
        check_nt = set([sequence[pos].upper() for pos in range(0, len(sequence))])
        if check_nt <= set(dna_base_options):
            print('Confirmed sequence contains exclusive DNA basepairs')
        else:
            raise ValueError('%s is not a valid nucleotide for a DNA object.' % check_nt)
        
        self.sequence = sequence.upper()
    
    
    def complement(self) -> str:
        '''
        input:
        description:
        return:
        '''
        dna_pairs = {'A':'T', 'T':'A', 'G':'C', 'C':'G'}

        complement_list = [dna_pairs[base] for base in list(self.sequence)]
        return ''.join(complement_list)
    
    def reverse(self) -> str:
        '''
        input:
        description:
        return:
        '''
        reverse_list = [base for base in list(self.sequence)[::-1]]
        return ''.join(reverse_list)
    
    def reverse_complement(self) -> str:
        '''
        input:
        description:
        return:
        '''
        dna_pairs = {'A':'T', 'T':'A', 'G':'C', 'C':'G'}

        reverse_list = [base for base in list(self.sequence)[::-1]]
        rc_list = [dna_pairs[base] for base in reverse_list]
        return ''.join(rc_list)
    
    def transcribe(self) -> Rna:
        '''
        input:
        description:
        return:
        '''
        return Rna(self.sequence.replace('T', 'U'))


    
    