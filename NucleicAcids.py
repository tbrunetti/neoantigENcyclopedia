from typing import *

class Rna:
    
    rna_pairs = {'A':'U', 'U':'A', 'G':'C', 'C':'G'}
    
    RNA_NUCLEOTIDES = ('A', 'U', 'G', 'C')
    
    def __init__(self, seq : str, rna_base_options = RNA_NUCLEOTIDES) -> None:
        check_nt = set([seq[pos].upper() for pos in range(0, len(seq))])
        if check_nt <= set(rna_base_options):
            print('Confirmed sequence contains exclusive RNA basepairs')
        else:
            raise ValueError('%s is not a valid nucleotide for a RNA object.' % check_nt)
        
        self.seq = seq.upper()
        
    
class Dna:
    dna_pairs = {'A':'T', 'T':'A', 'G':'C', 'C':'G'}
    
    DNA_NUCLEOTIDES = ('A', 'T', 'G', 'C')
    
    def __init__(self, seq : str, dna_base_options = DNA_NUCLEOTIDES) -> None:
        check_nt = set([seq[pos].upper() for pos in range(0, len(seq))])
        if check_nt <= set(dna_base_options):
            print('Confirmed sequence contains exclusive DNA basepairs')
        else:
            raise ValueError('%s is not a valid nucleotide for a DNA object.' % check_nt)
        
        self.seq = seq.upper()
    
    
    def complement(self) -> str:
        '''
        input:
        description:
        return:
        '''
        complement_list = [dna_pairs[base] for base in list(self.seq)]
        return ''.join(complement_list)
    
    def reverse(self) -> str:
        '''
        input:
        description:
        return:
        '''
        reverse_list = [base for base in list(self.seq)[::-1]]
        return ''.join(reverse_list)
    
    def reverse_complement(self) -> str:
        '''
        input:
        description:
        return:
        '''
        reverse_list = [base for base in list(self.seq)[::-1]]
        rc_list = [dna_pairs[base] for base in reverse_list]
        return ''.join(rc_list)
    
    def transcribe(self) -> Rna:
        '''
        input:
        description:
        return:
        '''
        return Rna(self.seq.replace('T', 'U'))


    
    