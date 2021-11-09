from Bio import SeqIO
import pandas
import argparse
import sys
import os
sys.path.append(os.getcwd())
from NucleicAcids import *
from Codon import *
from funcsForRefs import *
from typing import *


def calculate_orfs():
    '''
    input:
    description:
    return:
    '''
    pass


def combine_data(junctionFiles : list) -> pandas.DataFrame:
    '''
    input:
    description:
    return:
    '''
    
    # read in all SJ.tab files detected by 2-pass STAR and get the columns: chrom, intron start, intron end, strand, novelty
    junctionDFs = [pandas.read_table(inputFile, delim_whitespace = True, dtype = str, header = None, names = ['chrom', 'intron_start', 'intron_end', 'strandSTAR', 'motif', 'annotStatus'], 
                                   usecols = [0, 1, 2, 3, 4, 5]) for inputFile in junctionFiles]
    
    # merge all SJ.tab files after column extraction and only keep unique ones
    junctions = pandas.concat(junctionDFs, ignore_index = True)
    junctions.drop_duplicates(keep = 'last', inplace = True)
    
    '''
    rename columns values to be more meaningful:
        strand -> 0 = ud (undetermined), 1 = + strand, 2 = - strand
        motif -> 0 = nc (non-canonical), 1 = GT/AG, 2 = CT/AC, 3 = GC/AG, 4 = CT/GC, 5 = AT/AC, 6 = GT/AT
        annotStatus -> 0 = unannotated/novel, 1 = annotated
    '''
    junctions.replace({'strandSTAR': {'0':'ud' , '1': '+', '2':'-'}, 
                       'motif':{'0': 'nc', '1':'GT/AG', '2':'CT/AC', '3':'GC/AG', '4':'CT/GC', '5':'AT/AC', '6':'GT/AT'},
                       'annotStatus':{'0':'novel', '1':'annotated'}}, 
                       inplace = True)
    
    return junctions



def intron_retention(junctions : pandas.DataFrame, spladderOut : str, annotFilter : bool, diff_exp : str, pval_adj = float, outdir = str):
    '''
    input:
    description:
    return:
    '''
    
    # read in data from spladder and then rename columns for df merge
    as_events = pandas.read_csv(spladderOut, compression= 'gzip', sep = '\t', dtype = str)
    as_events.rename({'contig':'chrom'}, axis = 'columns', inplace = True)
    
    # merge data frame (ignore strand since some discrpenacies with undetermined in star vs determined in spladder)
    info_combine = as_events.merge(junctions, how = 'left', on = ['chrom', 'intron_start', 'intron_end'])
    # info_combine.to_csv("testing_introns.tsv", sep = "\t", index = False)
    
    # check if should filter based on annotation of junction being novel or annotated
    if annotFilter:
        print('Using junctions that are considered novel for peptide generation')
        filtered_junctions = info_combine.loc[info_combine['annotStatus'] == 'novel']
        filtered_junctions['event_id'] = filtered_junctions['event_id'].str.replace('.', '_')
    else:
        print('Using all annotations -- novel and annotated for peptide generation')
        filtered_junctions = info_combine
    
    # free up memory
    del info_combine
        
    # check if using differential ASE results
    if diff_exp != None:
        de_results = pandas.read_csv(diff_exp, sep = '\t')
        significant = de_results.loc[de_results['p_val_adj'] < 0.05]
        events_of_interest = significant.merge(filtered_junctions, on = 'event_id', how = 'inner')
        # TO DO: confirm gene and gene_name are the same
    else:
        events_of_interest = filtered_junctions
    
    # free up memory
    del filtered_junctions
    
    # TO DO: check strand
    # TO DO: check if strandSTAR and strand are concordant
    #       if not concordant, check if strandSTAR is ud
    #            if ud, then use strand
    #       if concordant then continue
    # TO DO: extract genomic sequence after getting all ORFs (3) and windows and convert to Dna()
    # TO DO: once extracted and once strand is determined then check strand:
    #       if strand is (-), then reverse compliment the extracted fasta sequence
    #       if strand is (+), nothing needs to be done
    
    
    
    # get psi info columns
    events_of_interest.filter(regex = '.psi', axis = 1)
    # get read coverage info columns
    events_of_interest.filter(regex = '_cov', axis = 1)
    # get mean event count and mean gene expression columns
    events_of_interest.filter(regex = '^mean_', axis = 1)
    # get log2FC of differential testing in event counts and gene expression levels columns
    events_of_interest.filter(regex = '^log2FC_', axis = 1)
   

def exon_skip():
    '''
    input:
    description:
    return:
    '''
    pass

def three_prime_alt():
    '''
    input:
    description:
    return:
    '''
    pass

def five_prime_alt():
    '''
    input:
    description:
    return:
    '''
    pass

def mutex():
    '''
    input:
    description:
    return:
    '''
    pass

def multiExon_skip():
    '''
    input:
    description:
    return:
    '''
    pass


if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(description= 'identification of peptides derived from alternative splice events from neojunctions')
    
    parser.add_argument('--SJfiles', required = True, action = 'extend', nargs= '+', type = str, help = 'list of all SJ.out.tab files from STAR two-pass mode')
    parser.add_argument('--deTestResults', type = str, default = None, help = 'spladder output of differentially expressed alternative splice events, if using results based on DE')
    parser.add_argument('--ase', type = str, help = 'file of type .txt.gz file of merged splice graphs output by spladder')
    parser.add_argument('--pvalAdj', type = float, default = 0.05, help = 'maximum [non-inclusive] pvalue adjusted cutoff to filter de results; only applicable if using --deTestResults')
    parser.add_argument('--outdir', type = str, default = os.getcwd(), help = 'directory to output results')
    parser.add_argument('--novel', action = 'store_true', type = bool, help = 'when this flag is specifided, it means to only return results where the splice junction is novel/unannotated \
                        otherwise, the default action is to return all alternative splice events')
    parser.add_argument('--fasta', type = str, help = 'fasta file of the version of the genome used to align and detect splice variants')
    args = parser.parse_args()
    
    codon_library = generate_codon_reference()    
    dna_fasta = SeqIO.index(args.fasta, 'fasta')
    junctions = combine_data(junctionFiles = args.SJfiles)
    intron_retention(junctions = junctions, spladderOut = args.ase, annotFilter = args.novel, diff_exp = args.deTestResults, pval_adj = args.pvalAdj, outdir = args.outdir)