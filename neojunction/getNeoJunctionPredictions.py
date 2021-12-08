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
import gc


def filtering_criteria(event_type :str, base_df : pandas.DataFrame, annot_filter : bool, diff_exp : Union[str, None], pval_adj : float)  -> pandas.DataFrame:
    '''
    input:
    description:
    return:
    '''
    
    # check if using only novel junctions
    if event_type == 'intron_retention':
            # check if should filter based on annotation of junction being novel or annotated
        if annot_filter:
            print('Using junctions that are considered novel for peptide generation')
            filtered_junctions = base_df.loc[base_df['annotStatus'] == 'novel']
            filtered_junctions['event_id'] = filtered_junctions['event_id'].str.replace('.', '_')
        else:
            print('Using all annotations -- novel and annotated for peptide generation')
            filtered_junctions = base_df
    
   
    elif event_type == 'exon_skip':
        # check if should filter based on annotation of junction being novel or annotated
        if annot_filter:
            print('Using junctions that are considered novel for peptide generation')
            filtered_junctions = base_df.loc[(base_df['annotStatus_pre'] == 'novel') | (base_df['annotStatus_aft'] == 'novel')]
            filtered_junctions['event_id'] = filtered_junctions['event_id'].str.replace('.', '_')
        else:
            print('Using all annotations -- novel and annotated for peptide generation')
            filtered_junctions = base_df
    
    elif event_type == 'three_prime_alt':
        # check if should filter based on annotation of junction being novel or annotated
        if annot_filter:
            print('Using junctions that are considered novel for peptide generation')
            filtered_junctions = base_df.loc[(base_df['annotStatus_alt1'] == 'novel') | (base_df['annotStatus_alt2'] == 'novel')]
            filtered_junctions['event_id'] = filtered_junctions['event_id'].str.replace('.', '_')
        else:
            print('Using all annotations -- novel and annotated for peptide generation')
            filtered_junctions = base_df
    


    # check if using differential ASE results
    if diff_exp != None:
        de_results = pandas.read_csv(diff_exp, sep = '\t')
        significant = de_results.loc[de_results['p_val_adj'] < pval_adj]
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
    

    return events_of_interest

@profile
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
        print("kmer is out of range for the full continuous region")
    
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
        print(rna.translate(codon_library, False))
        orf_ids['orf_{}_region'.format(orf)] = rna.translate(codon_library, False)
    
    peptide_bank[event_id] = orf_ids
    
    return peptide_bank
    
    

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



def intron_retention(junctions : pandas.DataFrame, spladderOut : str, outdir = str) -> None:
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
    
    events_of_interest = filtering_criteria(event_type = 'intron_retention', base_df = info_combine, annot_filter = args.novel, diff_exp = args.deTestResults, pval_adj = args.pvalAdj)
    
    # free up memory
    del info_combine
    gc.collect()
    
    # check if strandSTAR and strand are concordant
    #       if not concordant, check if strandSTAR is ud
    #            if ud, then use strand
    #       if concordant then continue
    for idx,row in events_of_interest.iterrows():
        if ((row['strandSTAR'] != row['strand']) & (row['strandSTAR'] == 'ud')):
            orfs = calculate_orfs(fasta = dna_fasta, kmer_length = args.kmer, strand = row['strand'], chrom = row['chrom'], flank_left_start = int(row['exon1_start']), flank_left_end = int(row['exon1_end']), flank_right_start = int(row['exon2_start']), flank_right_end = int(row['exon2_end']), ase_start = int(row['intron_start']), ase_end = int(row['intron_end']))
            peptides = translate_orfs(event_id = row['event_id'], orfs =  orfs, peptide_bank = peptides)
            
        elif ((row['strandSTAR'] != row['strand']) & (row['strandSTAR'] != 'ud')):
            print('star strand is {} and spladder strands is {}.  Skipping event.'.format(row['strandSTAR'] + row['strand']))
        
        elif (row['strandSTAR'] == row['strand']):
            orfs = calculate_orfs(fasta = dna_fasta, kmer_length = args.kmer, strand = row['strandSTAR'], chrom = row['chrom'], flank_left_start = int(row['exon1_start']), flank_left_end = int(row['exon1_end']), flank_right_start = int(row['exon2_start']), flank_right_end = int(row['exon2_end']), ase_start = int(row['intron_start']), ase_end = int(row['intron_end']))
            peptides = translate_orfs(event_id = row['event_id'], orfs = orfs, peptide_bank = peptides)
        
        del orfs
        gc.collect()

    # TO DO: calcuate percent overlap of kmer with flanking region
    
    
    
    # get psi info columns
    events_of_interest.filter(regex = '.psi', axis = 1)
    # get read coverage info columns
    events_of_interest.filter(regex = '_cov', axis = 1)
    # get mean event count and mean gene expression columns
    events_of_interest.filter(regex = '^mean_', axis = 1)
    # get log2FC of differential testing in event counts and gene expression levels columns
    events_of_interest.filter(regex = '^log2FC_', axis = 1)
   

def exon_skip(junctions : pandas.DataFrame, spladderOut : str, outdir = str) -> None:
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
    
    # generate columns to specify junctions between exon skips since star records intron junction not exon
    as_events['pre_intron_start'] = as_events['exon_pre_end'].astype(int) + 1
    as_events['pre_intron_end'] = as_events['exon_start'].astype(int) - 1
    as_events['aft_intron_start'] = as_events['exon_end'].astype(int) + 1
    as_events['aft_intron_end'] = as_events['exon_aft_start'].astype(int) - 1
    
    # merge data frame (ignore strand since some discrpenacies with undetermined in star vs determined in spladder)
    as_events['pre_intron_start'] = as_events['pre_intron_start'].astype(str)
    as_events['pre_intron_end'] = as_events['pre_intron_end'].astype(str)
    as_events['aft_intron_start'] = as_events['aft_intron_start'].astype(str)
    as_events['aft_intron_end'] = as_events['aft_intron_end'].astype(str)
    
    # rename to match the first intron start and end -- location of splice donor and acceptor sites
    tmp = junctions.rename({'intron_start':'pre_intron_start', 'intron_end':'pre_intron_end', 'strandSTAR':'strandSTAR_pre', 'motif':'motif_pre', 'annotStatus':'annotStatus_pre' }, axis = 'columns')
    info_combine = as_events.merge(tmp, how = 'left', on = ['chrom', 'pre_intron_start', 'pre_intron_end'])
    # rename to match the second intron start and end 
    tmp = junctions.rename({'intron_start':'aft_intron_start', 'intron_end':'aft_intron_end', 'strandSTAR':'strandSTAR_aft', 'motif':'mofif_aft', 'annotStatus':'annotStatus_aft'}, axis = 'columns')
    final = info_combine.merge(tmp, how = 'left', on = ['chrom', 'aft_intron_start', 'aft_intron_end'])

    # apply data filtering criteria
    events_of_interest = filtering_criteria(event_type = 'exon_skip', base_df = final, annot_filter = args.novel, diff_exp = args.deTestResults, pval_adj = args.pvalAdj)

    # free up memory
    del tmp, final, info_combine
    gc.collect()
    
    # check if strandSTAR and strand are concordant
    #       if not concordant, check if strandSTAR is ud
    #            if ud, then use strand
    #       if concordant then continue
    for idx,row in events_of_interest.iterrows():
        
        try:
            assert row['strandSTAR_pre'] == row['strandSTAR_aft'], 'There is a strand conflict at the following event_id: {}.  Skipping this event.'.format(row['event_id'])
        except  AssertionError as err:
            print(err)
            continue
            
        
        if ((row['strandSTAR_pre'] != row['strand']) & (row['strandSTAR_pre'] == 'ud')): # note it does not matter if we use strandSTAR_pre for strandStar_aft since we already checked the assertion that the STAR strands are concordant with each other
            orfs = calculate_orfs(fasta = dna_fasta, kmer_length = args.kmer, strand = row['strand'], chrom = row['chrom'], flank_left_start = int(row['exon_pre_start']), flank_left_end = int(row['exon_pre_end']), flank_right_start = int(row['exon_aft_start']), flank_right_end = int(row['exon_aft_end']), ase_start = int(row['exon_start']), ase_end = int(row['exon_end']))
            peptides = translate_orfs(event_id = row['event_id'], orfs =  orfs, peptide_bank = peptides)
        elif ((row['strandSTAR_pre'] != row['strand']) & (row['strandSTAR_pre'] != 'ud')):
            print('star strand is {} and spladder strands is {}.  Skipping event.'.format(row['strandSTAR_pre'] + row['strand']))
        
        
        elif (row['strandSTAR_pre'] == row['strand']):
            orfs = calculate_orfs(fasta = dna_fasta, kmer_length = args.kmer, strand = row['strandSTAR'], chrom = row['chrom'], flank_left_start = int(row['exon_pre_start']), flank_left_end = int(row['exon_pre_end']), flank_right_start = int(row['exon_aft_start']), flank_right_end = int(row['exon_aft_end']), ase_start = int(row['exon_start']), ase_end = int(row['exon_end']))
            peptides = translate_orfs(event_id = row['event_id'], orfs = orfs, peptide_bank = peptides)
        
        del orfs
        gc.collect() 
   
    # TO DO: calcuate percent overlap of kmer with flanking region
    
    
    # get psi info columns
    events_of_interest.filter(regex = '.psi', axis = 1)
    # get read coverage info columns
    events_of_interest.filter(regex = '_cov', axis = 1)
    # get mean event count and mean gene expression columns
    events_of_interest.filter(regex = '^mean_', axis = 1)
    # get log2FC of differential testing in event counts and gene expression levels columns
    events_of_interest.filter(regex = '^log2FC_', axis = 1)
   


@profile
def three_prime_alt(junctions : pandas.DataFrame, spladderOut : str, outdir = str) -> None:
    '''
    input:
    description:
    return:
    '''
    
    # dictionary that stores all ORF peptides for each event
    peptides_alt1 = {}
    peptides_alt2 = {}
    
    # read in data from spladder and then rename columns for df merge
    as_events = pandas.read_csv(spladderOut, compression= 'gzip', sep = '\t', dtype = str)
    as_events.rename({'contig':'chrom'}, axis = 'columns', inplace = True)
    
    # need to handle strands differently since this is specific to 3' alternative splice site    
    fwd_strand_events = as_events.loc[as_events['strand'] == '+']    
    rev_strand_events = as_events.loc[as_events['strand'] == '-']
    

    fwd_strand_events['alt1_intron_start'] = fwd_strand_events['exon_const_end'].astype(int) + 1
    fwd_strand_events['alt1_intron_end'] = fwd_strand_events['exon_alt1_start'].astype(int) - 1
    fwd_strand_events['alt2_intron_start'] = fwd_strand_events['exon_const_end'].astype(int) + 1
    fwd_strand_events['alt2_intron_end'] = fwd_strand_events['exon_alt2_start'].astype(int) - 1

 
    rev_strand_events['alt1_intron_start'] = rev_strand_events['exon_const_start'].astype(int) - 1
    rev_strand_events['alt1_intron_end'] = rev_strand_events['exon_alt1_end'].astype(int) + 1
    rev_strand_events['alt2_intron_start'] = rev_strand_events['exon_const_start'].astype(int) - 1
    rev_strand_events['alt2_intron_end'] = rev_strand_events['exon_alt2_end'].astype(int) + 1
    
    
    #dictionary specifity which columns require convert to type string to match type of star input file
    convert_types = {'alt1_intron_start': str,
                     'alt1_intron_end' : str,
                     'alt2_intron_start' : str,
                     'alt2_intron_end' : str
                     }
    
    
    fwd_strand_events = fwd_strand_events.astype(convert_types)
    rev_strand_events = rev_strand_events.astype(convert_types)

    # merge star intron annotation into fowrad strand events
    fwd_tmp = junctions.rename({'intron_start':'alt1_intron_start', 'intron_end':'alt1_intron_end', 'strandSTAR':'strandSTAR_alt1', 'motif':'motif_alt1', 'annotStatus':'annotStatus_alt1' }, axis = 'columns')
    fwd_alt1_info_combine = fwd_strand_events.merge(fwd_tmp, how = 'left', on = ['chrom', 'alt1_intron_start', 'alt1_intron_end'])
    
    fwd_tmp = junctions.rename({'intron_start':'alt2_intron_start', 'intron_end':'alt2_intron_end', 'strandSTAR':'strandSTAR_alt2', 'motif':'motif_alt2', 'annotStatus':'annotStatus_alt2' }, axis = 'columns')
    fwd_final = fwd_alt1_info_combine.merge(fwd_tmp, how = 'left', on = ['chrom', 'alt2_intron_start', 'alt2_intron_end'])

      
    # merge star intron annotation into reverse strand events
    rev_tmp = junctions.rename({'intron_start':'alt1_intron_end', 'intron_end':'alt1_intron_start', 'strandSTAR':'strandSTAR_alt1', 'motif':'motif_alt1', 'annotStatus':'annotStatus_alt1' }, axis = 'columns')
    rev_alt1_info_combine = rev_strand_events.merge(rev_tmp, how = 'left', on = ['chrom', 'alt1_intron_start', 'alt1_intron_end'])

    rev_tmp = junctions.rename({'intron_start':'alt2_intron_end', 'intron_end':'alt2_intron_start', 'strandSTAR':'strandSTAR_alt2', 'motif':'motif_alt2', 'annotStatus':'annotStatus_alt2' }, axis = 'columns')
    rev_final = rev_alt1_info_combine.merge(rev_tmp, how = 'left', on = ['chrom', 'alt2_intron_start', 'alt2_intron_end'])
    

    # merge dataframes back together
    combined_strand_events = pandas.concat([fwd_final, rev_final], axis = 0, ignore_index= True)
    
    # free memory
    del fwd_tmp, fwd_alt1_info_combine, rev_tmp, rev_alt1_info_combine, fwd_final, rev_final
    gc.collect()
    
    # apply data filtering criteria
    events_of_interest = filtering_criteria(event_type = 'three_prime_alt', base_df = combined_strand_events, annot_filter = args.novel, diff_exp = args.deTestResults, pval_adj = args.pvalAdj)

    # check if strandSTAR and strand are concordant
    #       if alt1 and alt2 are not concordant within just star, then skip\
    #       if not concordant between star and stran, check if strandSTAR is ud
    #            if ud, then use strand
    #       if concordant then continue
    
    for idx,row in events_of_interest.iterrows():
        
        try:
            assert row['strandSTAR_alt1'] == row['strandSTAR_alt2'], 'There is a strand conflict at the following event_id: {}.  Skipping this event.'.format(row['event_id'])
        except  AssertionError as err:
            print(err)
            continue
            
        # check strand for alt events since genomic positions are in reverse order for this alternative splicing event type
        
        if ((row['strandSTAR_alt1'] != row['strand']) & (row['strandSTAR_alt1'] == 'ud')): # note it does not matter if we use strandSTAR_pre for strandStar_aft since we already checked the assertion that the STAR strands are concordant with each other
            if row['strand'] == '+':
                orfs_alt1 = calculate_orfs(fasta = dna_fasta, kmer_length = args.kmer, strand = row['strand'], chrom = row['chrom'], flank_left_start = int(row['exon_const_start']), flank_left_end = int(row['exon_const_end']), flank_right_start = int(row['exon_alt1_start']), flank_right_end = int(row['exon_alt1_end']), ase_start = int(row['alt1_intron_start']), ase_end = int(row['alt1_intron_end']))
                peptides_alt1 = translate_orfs(event_id = row['event_id'], orfs =  orfs_alt1, peptide_bank = peptides_alt1)
                del orfs_alt1
                gc.collect()
                orfs_alt2 = calculate_orfs(fasta = dna_fasta, kmer_length = args.kmer, strand = row['strand'], chrom = row['chrom'], flank_left_start = int(row['exon_const_start']), flank_left_end = int(row['exon_const_end']), flank_right_start = int(row['exon_alt2_start']), flank_right_end = int(row['exon_alt2_end']), ase_start = int(row['alt2_intron_start']), ase_end = int(row['alt2_intron_end']))
                peptides_alt2 = translate_orfs(event_id = row['event_id'], orfs =  orfs_alt2, peptide_bank = peptides_alt2)
                del orfs_alt2
                gc.collect()
                
                
            elif row['strand'] =='-':
                orfs_alt1 = calculate_orfs(fasta = dna_fasta, kmer_length = args.kmer, strand = row['strand'], chrom = row['chrom'], flank_left_start = int(row['exon_alt1_end']), flank_left_end = int(row['exon_alt1_start']), flank_right_start = int(row['exon_const_start']), flank_right_end = int(row['exon_const_end']), ase_start = int(row['alt1_intron_end']), ase_end = int(row['alt1_intron_start']))
                peptides_alt1 = translate_orfs(event_id = row['event_id'], orfs =  orfs_alt1, peptide_bank = peptides_alt1)
                del orfs_alt1
                gc.collect()
                orfs_alt2 = calculate_orfs(fasta = dna_fasta, kmer_length = args.kmer, strand = row['strand'], chrom = row['chrom'], flank_left_start = int(row['exon_alt2_end']), flank_left_end = int(row['exon_alt2_start']), flank_right_start = int(row['exon_const_start']), flank_right_end = int(row['exon_const_end']), ase_start = int(row['alt2_intron_end']), ase_end = int(row['alt2_intron_start']))
                peptides_alt2 = translate_orfs(event_id = row['event_id'], orfs =  orfs_alt2, peptide_bank = peptides_alt2)
                del orfs_alt2
                gc.collect()


        elif ((row['strandSTAR_alt1'] != row['strand']) & (row['strandSTAR_alt1'] != 'ud')):
            print('star strand is {} and spladder strands is {}.  Skipping event.'.format(row['strandSTAR_alt1'] + row['strand']))
        
        
        elif (row['strandSTAR_alt1'] == row['strand']):
            
            if row['strand'] == '+':
                orfs_alt1 = calculate_orfs(fasta = dna_fasta, kmer_length = args.kmer, strand = row['strandSTAR_alt1'], chrom = row['chrom'], flank_left_start = int(row['exon_const_start']), flank_left_end = int(row['exon_const_end']), flank_right_start = int(row['exon_alt1_start']), flank_right_end = int(row['exon_alt1_end']), ase_start = int(row['alt1_intron_start']), ase_end = int(row['alt1_intron_end']))
                peptides_alt1 = translate_orfs(event_id = row['event_id'], orfs =  orfs_alt1, peptide_bank = peptides_alt1)
                del orfs_alt1
                gc.collect()
                orfs_alt2 = calculate_orfs(fasta = dna_fasta, kmer_length = args.kmer, strand = row['strandSTAR_alt2'], chrom = row['chrom'], flank_left_start = int(row['exon_const_start']), flank_left_end = int(row['exon_const_end']), flank_right_start = int(row['exon_alt2_start']), flank_right_end = int(row['exon_alt2_end']), ase_start = int(row['alt2_intron_start']), ase_end = int(row['alt2_intron_end']))
                peptides_alt2 = translate_orfs(event_id = row['event_id'], orfs =  orfs_alt2, peptide_bank = peptides_alt2)
                del orfs_alt2
                gc.collect()

            elif row['strand'] =='-':
                orfs_alt1 = calculate_orfs(fasta = dna_fasta, kmer_length = args.kmer, strand = row['strandSTAR_alt1'], chrom = row['chrom'], flank_left_start = int(row['exon_alt1_end']), flank_left_end = int(row['exon_alt1_start']), flank_right_start = int(row['exon_const_start']), flank_right_end = int(row['exon_const_end']), ase_start = int(row['alt1_intron_end']), ase_end = int(row['alt1_intron_start']))
                peptides_alt1 = translate_orfs(event_id = row['event_id'], orfs =  orfs_alt1, peptide_bank = peptides_alt1)
                del orfs_alt1
                gc.collect()
                orfs_alt2 = calculate_orfs(fasta = dna_fasta, kmer_length = args.kmer, strand = row['strandSTAR_alt2'], chrom = row['chrom'], flank_left_start = int(row['exon_alt2_end']), flank_left_end = int(row['exon_alt2_start']), flank_right_start = int(row['exon_const_start']), flank_right_end = int(row['exon_const_end']), ase_start = int(row['alt2_intron_end']), ase_end = int(row['alt2_intron_start']))
                peptides_alt2 = translate_orfs(event_id = row['event_id'], orfs =  orfs_alt2, peptide_bank = peptides_alt2)
                del orfs_alt2
                gc.collect()
   
    # TO DO: calcuate percent overlap of kmer with flanking region
    
    
    # get psi info columns
    events_of_interest.filter(regex = '.psi', axis = 1)
    # get read coverage info columns
    events_of_interest.filter(regex = '_cov', axis = 1)
    # get mean event count and mean gene expression columns
    events_of_interest.filter(regex = '^mean_', axis = 1)
    # get log2FC of differential testing in event counts and gene expression levels columns
    events_of_interest.filter(regex = '^log2FC_', axis = 1)
       
    
    
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
    parser.add_argument('--novel', action = 'store_true', help = 'when this flag is specifided, it means to only return results where the splice junction is novel/unannotated \
                        otherwise, the default action is to return all alternative splice events')
    parser.add_argument('--fasta', type = str, help = 'fasta file of the version of the genome used to align and detect splice variants')
    parser.add_argument('--kmer', type = int, help = 'kmers legnth to build peptide windows')
    parser.add_argument('--eventType', type = str, help= 'type of alternative splice event to analyze', choices = ['intron_retention', 'exon_skip', 'three_prime_alt', 'five_prime_alt', 'mutex'])
    parser.add_argument('--modules', default = os.path.join(os.getcwd(), '..', 'modules'), help = "full path to the modules directory location")
    args = parser.parse_args()
    
    codon_library = generate_codon_reference()    
    dna_fasta = SeqIO.index(args.fasta, 'fasta')
    junctions = combine_data(junctionFiles = args.SJfiles)
    
    if args.eventType == 'intron_retention':
        intron_retention(junctions = junctions, spladderOut = args.ase, outdir = args.outdir)
        
    if args.eventType == 'exon_skip':
        exon_skip(junctions = junctions, spladderOut = args.ase, outdir = args.outdir)

    if args.eventType == 'three_prime_alt':
        three_prime_alt(junctions = junctions, spladderOut = args.ase, outdir = args.outdir)