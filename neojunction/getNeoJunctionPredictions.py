import Bio
from Bio import SeqIO
import pandas
import argparse
import sys
import os
sys.path.append(os.path.join(os.getcwd(), '..', 'modules'))
from NucleicAcids import *
from Codon import *
from funcsForRefs import *
from typing import *


def calculate_orfs(fasta : Bio.File._IndexedSeqFileDict, kmer_length : int, strand : str, chrom : str, flank_left_start: int, flank_left_end: int, flank_right_start: int, flank_right_end : int, ase_start : int, ase_end : int) -> list[Dna]:
    '''
    input:
    description:
    return:
    '''
    
    def extract_continuity_regions() -> Tuple[list[int], Dna]:
        # create a new DNA that extracts and merges the continuous ORF regions into a new DNA sequence from genomic fasta
        left = Dna(fasta[chrom].seq[flank_left_start - 1 : flank_left_end + 1])
        event = Dna(fasta[chrom].seq[ase_start - 1 : ase_end + 1])
        right = Dna(fasta[chrom].seq[flank_right_start - 1 : flank_right_end + 1])

        # index where the start of the event and end of the event  occurs
        event_index = [len(left.sequence), len(left.sequence + event.sequence)-1]

        return  event_index, Dna(left.sequence + event.sequence + right.sequence)
   
   
    # TO DO: Find full length sequence to get all kmers and all ORFs
    orf_1_start = ((kmer_length - 1) * 3)  # number of bases to append upstream if 1st base of ase = first base in codon
    orf_2_start = ((kmer_length - 1) * 3) - 1  # number of bases to append upstream if 1st base of ase = second base in codon
    orf_3_start = ((kmer_length - 1) * 3) - 2 # number of bases to append upstream if 1st base of ase = last base in codon
    
    # confirmed these values
    # get orf ending based on if ase region is a multiple of 3 bases or not to append to last index of event
    if ((ase_end - ase_start) + 1) % 3 == 0:
        print('mod_result_0')
        orf_1_end = ((kmer_length - 1) * 3) 
        orf_2_end = ((kmer_length - 1) * 3) + 1
        orf_3_end = ((kmer_length - 1) * 3) + 2
    elif ((ase_end - ase_start) + 1) % 3 == 1:
        print('mod_result_1')
        orf_1_end = ((kmer_length - 1) * 3) + 2
        orf_2_end = ((kmer_length - 1) * 3) 
        orf_3_end = ((kmer_length - 1) * 3) + 1
    elif ((ase_end - ase_start) + 1) % 3 == 2:
        print('mod_result_2')
        orf_1_end = ((kmer_length - 1) * 3) + 1
        orf_2_end = ((kmer_length - 1) * 3) + 2
        orf_3_end = ((kmer_length - 1) * 3)         
    
    
   
    # TO DO: extract genomic sequence after getting all ORFs (3) and windows and convert to Dna() -- python Seq is 0-indexed, genomic coords are 1-indexed
    
    event_index, dna_continuous_seq = extract_continuity_regions()
    try:
        orf_1_region = Dna(dna_continuous_seq.sequence[event_index[0] - orf_1_start : event_index[1] + orf_1_end])
        orf_2_region = Dna(dna_continuous_seq.sequence[event_index[0] - orf_2_start : event_index[1] + orf_2_end])
        orf_3_region = Dna(dna_continuous_seq.sequence[event_index[0] - orf_3_start : event_index[1] + orf_3_end])
    except IndexError:
        print("kmer is out of range for the full continous region")
    
    #orf_1_region = Dna(fasta[chrom].sequence[orf_1_start - 1 : orf_1_end + 1])
    #orf_2_region = Dna(fasta[chrom].sequence[orf_2_start - 1 : orf_2_end + 1])
    #orf_3_region = Dna(fasta[chrom].sequence[orf_3_start - 1 : orf_3_end + 1])
    
    # TO DO: once extracted and once strand is determined then check strand:
    # if strand is (-), then reverse complement the extracted fasta sequence
    # if strand is (+), nothing needs to be done
    if strand == '-':
        orf_1_region = Dna(orf_1_region.reverse_complement())
        orf_2_region = Dna(orf_2_region.reverse_complement())
        orf_3_region = Dna(orf_3_region.reverse_complement())
    
    return [orf_1_region, orf_2_region, orf_3_region]
        


def translate_orfs(event_id: str, orfs : list[Dna], peptide_bank : Dict[str, Dict[str,str]]) -> Dict[str, Dict[str, str]]:
    rna_list = [] # a list of Rna objects
    orf_ids = {}
    for dna in orfs:
        rna_list.append(dna.transcribe()) # transcribe() returns an Rna object
        
    for orf, rna in enumerate(rna_list):
        print((orf, rna))
        orf_ids['orf_{}_region'.format(orf)] = rna.sequence.translate(codon_library)
    
    return peptide_bank[event_id : orf_ids]
    
    

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



def intron_retention(junctions : pandas.DataFrame, spladderOut : str, annotFilter : bool, diff_exp : str, pval_adj = float, outdir = str) -> None:
    '''
    input:
    description:
    return:
    '''
    
    # dictionary that stores all ORF peptides for each event
    peptides = {}
    
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
        # confirm gene names and ID match between two dataframes
        if events_of_interest['gene'].equals(events_of_interest['gene_name']):
            print('CONFIRMED -- all gene names match')
        else:
            # if gene and gene_name are not the same error is raised with list of mismatches
            mismatches = list(events_of_interest['event_id'].loc[events_of_interest['gene'] != events_of_interest['gene_name']])
            raise ValueError('FAILED! Dataframe merge has gene mismatches between spladder input and STAR junctions. \
                             Check the following event ids: {}'.format(','.join(mismatches)))
    else:
        events_of_interest = filtered_junctions
    
    # free up memory
    del filtered_junctions
    
    # TO DO: check strand
    # TO DO: check if strandSTAR and strand are concordant
    #       if not concordant, check if strandSTAR is ud
    #            if ud, then use strand
    #       if concordant then continue
    for idx,row in events_of_interest.iterrows():
        if ((row['strandSTAR'] != row['strand']) & (row['strandSTAR'] == 'ud')):
            orfs = calculate_orfs(fasta = dna_fasta, kmer_length = 8, strand = row['strand'], chrom = row['chrom'], flank_left_start = row['exon1_start'], flank_left_end = row['exon1_end'], flank_right_start = row['exon2_start'], flank_right_end = row['exon2_end'], ase_start = row['intron_start'], ase_end = row['intron_end'])
            peptides = translate_orfs(event_id = row['event_id'], orfs =  orfs, peptide_bank = peptides)
        elif ((row['strandSTAR'] != row['strand']) & (row['strandSTAR'] != 'ud')):
            print('star strand is {} and spladder strands is {}.  Skipping event.'.format(row['strandSTAR'] + row['strand']))
        
        elif (row['strandSTAR'] == row['strand']):
            orfs = calculate_orfs(fasta = dna_fasta, kmer_length = 8, strand = row['strandSTAR'], chrom = row['chrom'], flank_left_start = row['exon1_start'], flank_left_end = row['exon1_end'], flank_right_start = row['exon2_start'], flank_right_end = row['exon2_end'], ase_start = row['intron_start'], ase_end = row['intron_end'])
            peptides = translate_orfs(event_id = row['event_id'], orfs = orfs, peptide_bank = peptides)

    # TO DO: calcuate percent overlap of kmer with flanking region
    
    
    
    # get psi info columns
    events_of_interest.filter(regex = '.psi', axis = 1)
    # get read coverage info columns
    events_of_interest.filter(regex = '_cov', axis = 1)
    # get mean event count and mean gene expression columns
    events_of_interest.filter(regex = '^mean_', axis = 1)
    # get log2FC of differential testing in event counts and gene expression levels columns
    events_of_interest.filter(regex = '^log2FC_', axis = 1)
   

def exon_skip(junctions : pandas.DataFrame, spladderOut : str, annotFilter : bool, diff_exp : str, pval_adj = float, outdir = str) -> None:
    '''
    input:
    description:
    return:
    '''
    # dictionary that stores all ORF peptides for each event
    peptides = {}
    
    # read in data from spladder and then rename columns for df merge
    as_events = pandas.read_csv(spladderOut, compression= 'gzip', sep = '\t', dtype = str)
    as_events.rename({'contig':'chrom'}, axis = 'columns', inplace = True)

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
    parser.add_argument('--kmer', type = int, help = 'kmers legnth to build peptide windows')
    parser.add_argument('--eventType', type = str, help= 'type of alternative splice event to analyze', choices = ['intron_retention', 'exon_skip', 'three_prime_alt', 'five_prime_alt', 'mutex'])
    args = parser.parse_args()
    
    codon_library = generate_codon_reference()    
    dna_fasta = SeqIO.index(args.fasta, 'fasta')
    junctions = combine_data(junctionFiles = args.SJfiles)
    
    if args.eventType == 'intron_retention':
        intron_retention(junctions = junctions, spladderOut = args.ase, annotFilter = args.novel, diff_exp = args.deTestResults, pval_adj = args.pvalAdj, outdir = args.outdir)