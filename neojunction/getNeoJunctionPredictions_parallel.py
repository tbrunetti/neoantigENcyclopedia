from multiprocessing import process
from subprocess import CalledProcessError
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
import math
import multiprocessing
import numpy
import time

def flags(final_df : pandas.DataFrame, kmer_length : int) -> pandas.DataFrame:
    '''
    input:
    description:
    return:
    '''
    final_df[['constant_region_premature_stop', 'short_ase_length', 'ase_stop_encountered']] = bool, bool, bool
    final_df['constant_region_premature_stop'] = final_df['upstream_region_translated_frame'].str.contains("\*")
    final_df['short_ase_length'] = final_df['translated_ase_peptide'].str.len() < kmer_length
    final_df['ase_stop_encountered'] = final_df['translated_ase_peptide'].str.contains("\*")
    
    strict_clean_flags_df = final_df.loc[(final_df['constant_region_premature_stop'] == False) & (final_df['short_ase_length'] == False) & (final_df['blast_expected_value'].astype(float) < 0.0001)]

    strict_clean_flags_df.to_csv('final_dataframe_test_cleaned_flags.txt', sep = "\t")
    
    fasta_netmhc_file = open('intron_retentions_netmhc_file.fasta', 'w')
    
    for idx,row in strict_clean_flags_df.iterrows():
        fasta_netmhc_file.write('>{}\n{}\n'.format(idx, row['translated_ase_peptide']))
        fasta_netmhc_file.flush()
    
    fasta_netmhc_file.close()
    

    return final_df


@profile
def blast_search(event_id: str, blast : str, db : str, fasta : Bio.File._IndexedSeqFileDict, chrom : str, const_region_start : int, const_region_end : int, strand : str, gene : str, blast_tmp_file : int):
    import subprocess
    import re

    '''
    TO DO: check if file exists to make sure to warn user it will be overwritten in even they happen to have this file in their space
    '''
    
    if strand == '+':
        print(event_id, fasta[chrom].seq[const_region_start - 1 : const_region_end], const_region_start, const_region_end, chrom, const_region_start -1)
        try:
            reg_1 = Dna(fasta[chrom].seq[const_region_start - 1 : const_region_end]).transcribe().translate(codon_library, True)
        except ValueError:
            reg_1 = ''
        try:
            reg_2 = Dna(fasta[chrom].seq[const_region_start - 2 : const_region_end]).transcribe().translate(codon_library, True)
        except ValueError: 
            reg_2 = ''
        try:
            reg_3 = Dna(fasta[chrom].seq[const_region_start - 3 : const_region_end]).transcribe().translate(codon_library, True)
        except ValueError:
            reg_3 = ''

    elif strand == '-': # reverse complements region before translation
        try:
            reg_1 = Dna(Dna(fasta[chrom].seq[const_region_start - 1 : const_region_end]).reverse_complement()).transcribe().translate(codon_library, True)
        except ValueError:
            reg_1 = ''
        try:
            reg_2 = Dna(Dna(fasta[chrom].seq[const_region_start - 2 : const_region_end]).reverse_complement()).transcribe().translate(codon_library, True)
        except ValueError:
            reg_2 = ''
        try:
            reg_3 = Dna(Dna(fasta[chrom].seq[const_region_start - 3 : const_region_end]).reverse_complement()).transcribe().translate(codon_library, True)
        except ValueError:
            reg_3 = ''


    # generate an ORF
    blastFile = open(str(blast_tmp_file) + '_parallel_input_search.fasta', 'w')
    blastFile.write('>{}\n{}\n'.format(reg_1, reg_1))
    blastFile.write('>{}\n{}\n'.format(reg_2, reg_2))
    blastFile.write('>{}\n{}'.format(reg_3, reg_3))
    blastFile.flush() # flush the buffer
    blastFile.close() # no longer need to write to it, so close file

    command = '{} -query {} -db {} -subject_besthit -outfmt "6 delim=, std salltitles"'.format(blast, blastFile.name, db)
    try:
        blast_results = subprocess.check_output([command], shell = True) # returns a byte level blast result that should be decoded into a string
        print(blast_results)
    except CalledProcessError:
        print("Error in subrocess call!  Skipping result")
        return ('none found', 'none found', 'none found', 'none found')

    # check that return code is not anything other than 0
    '''
    try:
        blast_results.check_returncode
    except CalledProcessError as e:
        print('There was a runtime error in your blast search. \n {}'.format(e))
    '''
        
    # decode blast results and extract the matching protein of found frame, corresponding preotein ID, and gene symbol 
    #frame_matches = re.findall(r'(.+?),(.+?),.+?GN=(.+?)\s', blast_results.decode()) # returns a list of tuples each of length 3
    frame_matches = re.findall(r'(.+?),(.+?),.+?,.+?,.+?,.+?,.+?,.+?,.+?,.+?,(.+?),.+?GN=(.+?)\s', blast_results.decode()) # returns a list of tuples each of length 4
    print(frame_matches, gene, event_id)
    
    # find matching frame of gene name listed relative to upstream constant region identified
    all_correct_frames = [blastMatch for blastMatch in frame_matches if gene in blastMatch]
    correct_frame = list(set(all_correct_frames))
    print(correct_frame)
    
    if len(correct_frame) == 0: # means none of the frame match the gene identified in spladder df
        del blastFile
        gc.collect()
        return ('none found', 'none found', 'none found', 'none found')
    elif len(correct_frame) == 1: # means only a single frame matches the gene idenfied in spladder
        del blastFile
        gc.collect()
        return correct_frame[0]
    else: # means multiple frames matches the gene identified in spladder -- pick one with best e-value from BLAST
        del blastFile
        gc.collect()
        evals = []
        for results in correct_frame:
            evals.append(float(results[2]))
        return correct_frame[evals.index(min(evals))] 
    

def filtering_criteria(event_type :str, base_df : pandas.DataFrame, geneMatch : str, annot_filter : bool, diff_exp : Union[str, None], pval_adj : float)  -> pandas.DataFrame:
    '''
    input:
    description:
    return:
    '''
    
    def format_gtf_column():
        tmp_pair = {}
        with open(geneMatch, 'r') as nameInfo:
            for line in nameInfo:
                if (line.split(':')[0]) in tmp_pair:
                    tmp_pair[line.split(':')[0]].append(line.split(':')[1].replace('gene_id \"', '').replace('\";', '').replace('gene_name \"', '').strip())
                else:
                    tmp_pair[line.split(':')[0]] = [line.split(':')[1].replace('gene_id \"', '').replace('\";', '').replace('gene_name \"', '').strip()]

        geneMap = pandas.DataFrame(tmp_pair.values(), columns = ['gene_name', 'symbol'])
        #free memory
        del tmp_pair
        gc.collect()
        geneMap.drop_duplicates(keep = 'first', inplace = True)

        return geneMap
        
    
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
    
    elif (event_type == 'three_prime_alt') | (event_type == 'five_prime_alt'):
        # check if should filter based on annotation of junction being novel or annotated
        if annot_filter:
            print('Using junctions that are considered novel for peptide generation')
            filtered_junctions = base_df.loc[(base_df['annotStatus_alt1'] == 'novel') | (base_df['annotStatus_alt2'] == 'novel')]
            filtered_junctions['event_id'] = filtered_junctions['event_id'].str.replace('.', '_')
        else:
            print('Using all annotations -- novel and annotated for peptide generation')
            filtered_junctions = base_df
    
    elif event_type == 'mutex':
        if annot_filter:
            print('Using junctions that are considered novel for peptide generation')
            filtered_junctions = base_df.loc[(base_df['annotStatus_exon1_pre'] == 'novel') | (base_df['annotStatus_exon1_aft'] == 'novel') | (base_df['annotStatus_exon2_pre'] == 'novel') | (base_df['annotStatus_exon2_aft'] == 'novel')]
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
    

    # perform blast frame
    geneMap = ''
    with open(geneMatch, 'r') as testLine:
        # means needs to be formatted since it is in grep output format
        if len(testLine.readline().split('\t')) == 1:
            geneMap = format_gtf_column()
            annotated_events_of_interest = events_of_interest.merge(geneMap, how = 'left', on = 'gene_name')
        # means already formatted and ready to merge
        else:
            geneMap = pandas.read_csv(geneMatch, sep = '\t',  header = None, names = ['gene_name', 'symbol']) 
            annotated_events_of_interest = events_of_interest.merge(geneMap, how = 'left', on = 'gene_name')

    return annotated_events_of_interest


def calculate_orfs(event_type :str, fasta : Bio.File._IndexedSeqFileDict, kmer_length : int, strand : str, chrom : str, flank_left_start: int, flank_left_end: int, flank_right_start: int, flank_right_end : int, ase_start : Union[int, None], ase_end : Union[int, None]) -> list[Dna]:
    '''
    input:
    description:
    return:
    '''
    
    def extract_continuity_regions() -> Tuple[list[int], Dna]:
        
        # create a new DNA that extracts and merges the continuous ORF regions into a new DNA sequence from genomic fasta
        if ase_end == None: # sometimes there will not be an end and should be translated through the whole region -- such as alt 5' and alt 3' events
            try:
                left = Dna(fasta[chrom].seq[flank_left_start - 1 : flank_left_end])
            except ValueError:
                left = ''
            
            try:
                right = Dna(fasta[chrom].seq[flank_right_start - 1 : flank_right_end])
            except ValueError:
                right = ''
            
            event_index = [len(left.sequence), len(right.sequence)]
            
            #if args.eventType == 'three_prime_alt':            
            #    event_index = [len(left.sequence), 0] # len(left.sequence) would be the first base of the right sequence
            
            #elif args.eventType == 'five_prime_alt':
            #    event_index = [len(right.sequence), 0]
            
            return event_index, Dna(left.sequence + right.sequence)
       
        
        elif args.eventType == 'exon_skip':
            try:
                left = Dna(fasta[chrom].seq[flank_left_start - 1 : flank_left_end])
            except ValueError:
                left = ''
            
            try:
                event_left = Dna(fasta[chrom].seq[ase_start : ase_start + 3]) # always 3 nucleotides last full codon of exon is max
            except ValueError:
                event_left = ''
            
            try:
                event_right = Dna(fasta[chrom].seq[ase_end - 4 : ase_end - 1]) # always 3 nucleotides first full codon of exon is max
            except ValueError:
                event_right = ''
                
            try:
                right = Dna(fasta[chrom].seq[flank_right_start - 1 : flank_right_end])
            except ValueError:
                right = ''
                
            event_index = [len(left.sequence), len(left.sequence + event_left.sequence + event_right.sequence)-1]
            return event_index, Dna(left.sequence + event_left.sequence + event_right.sequence + right.sequence)

        
        else:
            try:
                left = Dna(fasta[chrom].seq[flank_left_start - 1 : flank_left_end])
            except ValueError:
                left = ''
            
            try:
                event = Dna(fasta[chrom].seq[ase_start - 1 : ase_end])
            except ValueError:
                event = ''
                
            try:
                right = Dna(fasta[chrom].seq[flank_right_start - 1 : flank_right_end])
            except:
                right = ''

        # index where the start of the event and end of the event  occurs
        event_index = [len(left.sequence), len(left.sequence + event.sequence)-1]

        return  event_index, Dna(left.sequence + event.sequence + right.sequence)
   
   
    # TO DO: Find full length sequence to get all kmers and all ORFs
    orf_1_start = ((kmer_length - 1) * 3) # number of bases to append upstream if 1st base of ase = first base in codon = mod0
    orf_2_start = ((kmer_length - 1) * 3) - 1  # number of bases to append upstream if 1st base of ase = second base in codon = mod2
    orf_3_start = ((kmer_length - 1) * 3) - 2 # number of bases to append upstream if 1st base of ase = last base in codon = mod1
    
    # confirmed these values
    # get orf ending based on if ase region is a multiple of 3 bases or not to append to last index of event
    if ase_end != None: # if this holds true, it means the ase is not in between two constant regions
        if ((ase_end - ase_start) + 1) % 3 == 0:
            print('mod_result_0')
            orf_1_end = ((kmer_length - 1) * 3) 
            orf_3_end = ((kmer_length - 1) * 3) + 2
            orf_2_end = ((kmer_length - 1) * 3) + 1

        elif ((ase_end - ase_start) + 1) % 3 == 1:
            print('mod_result_1')
            orf_3_end = ((kmer_length - 1) * 3) + 2
            orf_1_end = ((kmer_length - 1) * 3) 
            orf_2_end = ((kmer_length - 1) * 3) + 1
           
        elif ((ase_end - ase_start) + 1) % 3 == 2:
            print('mod_result_2')
            orf_2_end = ((kmer_length - 1) * 3) + 1
            orf_3_end = ((kmer_length - 1) * 3) + 2
            orf_1_end = ((kmer_length - 1) * 3)  
         

   
    # TO DO: extract genomic sequence after getting all ORFs (3) and windows and convert to Dna() -- python Seq is 0-indexed, genomic coords are 1-indexed
    event_index, dna_continuous_seq = extract_continuity_regions()
    
    try:
        # means that ase event is located between two constant regions (i.e. intron retention and mutex)
        if ase_end != None: 
            # if upstream and downstream exon is smaller than kmer window take full region (+ strand)
            if (event_index[0] < max(orf_1_start, orf_2_start, orf_3_start)) & ((len(dna_continuous_seq.sequence) - event_index[1]+1) < max(orf_1_end, orf_2_end, orf_3_end)) & (strand == '+'):
                orf_1_region = Dna(dna_continuous_seq.sequence)
                orf_2_region = Dna(dna_continuous_seq.sequence[1:])
                orf_3_region = Dna(dna_continuous_seq.sequence[2:])
                check_frame_region = math.floor(len(dna_continuous_seq.sequence[:event_index[0]])/3)
            # if upstream and downstream exon is smaller than kmer window take full region (- strand)
            elif (event_index[0] < max(orf_1_start, orf_2_start, orf_3_start)) & ((len(dna_continuous_seq.sequence) - event_index[1]+1) < max(orf_1_end, orf_2_end, orf_3_end)) & (strand == '-'):
                orf_1_region = Dna(dna_continuous_seq.sequence)
                orf_2_region = Dna(dna_continuous_seq.sequence[:-1])
                orf_3_region = Dna(dna_continuous_seq.sequence[:-2])
                check_frame_region = math.floor(len(dna_continuous_seq.sequence[event_index[1]+1:])/3)

            # if left exon is smaller than kmer widown but right flank exceeds kmer window (+ strand)
            elif (event_index[0] < max(orf_1_start, orf_2_start, orf_3_start)) & ((len(dna_continuous_seq.sequence) - event_index[1]+1) > max(orf_1_end, orf_2_end, orf_3_end)) & (strand == '+'):
                orf_1_region = Dna(dna_continuous_seq.sequence[:event_index[1] + orf_1_end])
                orf_2_region = Dna(dna_continuous_seq.sequence[1:event_index[1] + orf_2_end])
                orf_3_region = Dna(dna_continuous_seq.sequence[2:event_index[1] + orf_3_end]) 
                check_frame_region = math.floor(len(dna_continuous_seq.sequence[:event_index[0]])/3)
            # if left exon is smaller than kmer widown but right flank exceeds kmer window (- strand)
            elif (event_index[0] < max(orf_1_start, orf_2_start, orf_3_start)) & ((len(dna_continuous_seq.sequence) - event_index[1]+1) > max(orf_1_end, orf_2_end, orf_3_end)) & (strand == '-'):
                orf_1_region = Dna(dna_continuous_seq.sequence[:event_index[1] + orf_1_end])
                orf_2_region = Dna(dna_continuous_seq.sequence[1:event_index[1] + orf_2_end])
                orf_3_region = Dna(dna_continuous_seq.sequence[2:event_index[1] + orf_3_end])    
                check_frame_region = kmer_length - 1
            # if right exon is smaller than kmer widown but left flank exceeds kmer window (+ strand)
            elif (event_index[0] > max(orf_1_start, orf_2_start, orf_3_start)) & ((len(dna_continuous_seq.sequence) - event_index[1]+1) < max(orf_1_end, orf_2_end, orf_3_end)) & (strand == '+'):
                orf_1_region = Dna(dna_continuous_seq.sequence[event_index[0] - orf_1_start:])
                orf_2_region = Dna(dna_continuous_seq.sequence[event_index[0] - orf_2_start:])
                orf_3_region = Dna(dna_continuous_seq.sequence[event_index[0] - orf_3_start:])
                check_frame_region = kmer_length - 1
            # if right exon is smaller than kmer widown but left flank exceeds kmer window (- strand)
            elif (event_index[0] > max(orf_1_start, orf_2_start, orf_3_start)) & ((len(dna_continuous_seq.sequence) - event_index[1]+1) < max(orf_1_end, orf_2_end, orf_3_end)) & (strand == '-'):
                orf_1_region = Dna(dna_continuous_seq.sequence[event_index[0] - orf_1_start:])
                orf_2_region = Dna(dna_continuous_seq.sequence[event_index[0] - orf_2_start:-1])
                orf_3_region = Dna(dna_continuous_seq.sequence[event_index[0] - orf_3_start:-2])
                check_frame_region = math.floor(len(dna_continuous_seq.sequence[event_index[1]+1:])/3)
            # means both flanking exon exceed kmer windows
            else: 
                orf_1_region = Dna(dna_continuous_seq.sequence[event_index[0] - orf_1_start:event_index[1] + 1 + orf_1_end])
                orf_2_region = Dna(dna_continuous_seq.sequence[event_index[0] - orf_2_start:event_index[1] + 1 +  orf_2_end])
                orf_3_region = Dna(dna_continuous_seq.sequence[event_index[0] - orf_3_start:event_index[1] + 1 + orf_3_end])
                check_frame_region = kmer_length - 1   # number of amino acids to check upstream as a substring

        # means only the start or end is a constant regions and the remainder of the orf window should fully be translated through the end of the sequence (ex: alt 3', alt 5')
        else:
            if ((strand == '+') & (event_type == 'three_prime_alt')):
                # if that number of bases upstream from max start is less than the kmer window calculated start, then translate the whole thing and then shift the reading frame by 1 base for each orf
                if max(orf_1_start, orf_2_start, orf_3_start) > event_index[0]: 
                    orf_1_region = Dna(dna_continuous_seq.sequence)
                    orf_2_region = Dna(dna_continuous_seq.sequence[1:])
                    orf_3_region = Dna(dna_continuous_seq.sequence[2:])
                # otherwise, use the calculated start and go all the way through the end 
                else:
                    orf_1_region = Dna(dna_continuous_seq.sequence[event_index[0] - orf_1_start :])
                    orf_2_region = Dna(dna_continuous_seq.sequence[event_index[0] - orf_2_start :])
                    orf_3_region = Dna(dna_continuous_seq.sequence[event_index[0] - orf_3_start :])
            
            elif ((strand == '-') & (event_type == 'five_prime_alt')):
                 # TO DO check window
                if (max(orf_1_start, orf_2_start, orf_3_start) + event_index[1]) > len(dna_continuous_seq.sequence):
                    orf_1_region = Dna(dna_continuous_seq.sequence)
                    orf_2_region = Dna(dna_continuous_seq.sequence[:-1])
                    orf_3_region = Dna(dna_continuous_seq.sequence[:-2])
                # otherwise, use the calculated start and go all the way through the end 
                else:
                    orf_1_region = Dna(dna_continuous_seq.sequence[event_index[0] - orf_1_start :])
                    orf_2_region = Dna(dna_continuous_seq.sequence[event_index[0] - orf_2_start :-1])
                    orf_3_region = Dna(dna_continuous_seq.sequence[event_index[0] - orf_3_start :-2])
            
            elif ((strand == '-') & (event_type == 'three_prime_alt')):
                # TO DO check window
                if (max(orf_1_start, orf_2_start, orf_3_start) + event_index[0]) > len(dna_continuous_seq.sequence):
                    orf_1_region = Dna(dna_continuous_seq.sequence)
                    orf_2_region = Dna(dna_continuous_seq.sequence[:-1])
                    orf_3_region = Dna(dna_continuous_seq.sequence[:-2])
                # otherwise, use the calculated start and go all the way through the end 
                else:
                    orf_1_region = Dna(dna_continuous_seq.sequence[:event_index[0] + orf_1_start])
                    orf_2_region = Dna(dna_continuous_seq.sequence[1:event_index[0] + orf_2_start])
                    orf_3_region = Dna(dna_continuous_seq.sequence[2:event_index[0] + orf_3_start])
            elif ((strand == '+') & (event_type == 'five_prime_alt')):
                # if that number of bases upstream from max start is less than the kmer window calculated start, then translate the whole thing and then shift the reading frame by 1 base for each orf
                if max(orf_1_start, orf_2_start, orf_3_start) > event_index[1]: 
                    orf_1_region = Dna(dna_continuous_seq.sequence)
                    orf_2_region = Dna(dna_continuous_seq.sequence[1:])
                    orf_3_region = Dna(dna_continuous_seq.sequence[2:])
                # otherwise, use the calculated start and go all the way through the end 
                else:
                    orf_1_region = Dna(dna_continuous_seq.sequence[0:event_index[0] + orf_1_start])
                    orf_2_region = Dna(dna_continuous_seq.sequence[1:event_index[0] + orf_2_start+1])
                    orf_3_region = Dna(dna_continuous_seq.sequence[2:event_index[0] + orf_3_start+2])
                    
                print(dna_continuous_seq.sequence)
                print(orf_1_region.sequence, orf_2_region.sequence, orf_3_region.sequence)

    except IndexError:
        print("kmer is out of range for the full continuous region")
        
    # TO DO: once extracted and once strand is determined then check strand:
    # if strand is (-), then reverse complement the extracted fasta sequence
    # if strand is (+), nothing needs to be done
    if strand == '-':
        orf_1_region = Dna(orf_1_region.reverse_complement())
        orf_2_region = Dna(orf_2_region.reverse_complement())
        orf_3_region = Dna(orf_3_region.reverse_complement())
    
    return [orf_1_region, orf_2_region, orf_3_region], check_frame_region
        


def translate_orfs(event_id: str, orfs : list[Dna], peptide_bank : Dict[str, Dict[str,str]], upstream_const : str, correct_frame_only : bool, frame_check: int) -> Dict[str, Dict[str, str]]:
    rna_list = [] # a list of Rna objects
    orf_ids = {}
    matching_frames = {}
    for dna in orfs:
        rna_list.append(dna.transcribe()) # transcribe() returns an Rna object
        
    for orf, rna in enumerate(rna_list):
        print((orf, rna))
        print(event_id, rna.translate(codon_library, False))
        if correct_frame_only:
            # return the highest index where contsant upstream region was found
            matching_frames[orf] = upstream_const.rfind(rna.translate(codon_library, False)[:frame_check])
            orf_ids['orf_{}_region'.format(orf)] = rna.translate(codon_library, False)
        else:
            orf_ids['orf_{}_region'.format(orf)] = rna.translate(codon_library, False)
            
    if correct_frame_only:
        # confirm all keys are not -1
        try:
            # if not all -1, find the match closest to the ase event and select that reading frame
            assert max(matching_frames.values()) != -1
            correct_orf_key = max(matching_frames, key=matching_frames.get)
            peptide_bank[event_id] = orf_ids['orf_{}_region'.format(correct_orf_key)]
        except AssertionError:
            # if all are -1, then reduce frame match size by 1 and try again
            frame_check = frame_check - 1
            for orf, rna in enumerate(rna_list):
                matching_frames[orf] = upstream_const.rfind(rna.translate(codon_library, False)[:frame_check])
                orf_ids['orf_{}_region'.format(orf)] = rna.translate(codon_library, False)
                
            try: 
                assert max(matching_frames.values()) != -1
                correct_orf_key = max(matching_frames, key=matching_frames.get)
                peptide_bank[event_id] = orf_ids['orf_{}_region'.format(correct_orf_key)]
            # if after decreasing the constant region frame length by 1 and there are still no matches, then return that peptide does not match any reading frames
            except AssertionError:
                print('For event id, {}, there are no valid matching open reading frames with know protein accession ids.'.format(event_id))
                peptide_bank[event_id] = 'No open reading frame detected'
    
    else:
        peptide_bank[event_id] = orf_ids
    
    
    return peptide_bank
    
    

def combine_data(junctionFiles : list, annotation : str) -> pandas.DataFrame:
    '''
    input:
    description:
    return:
    '''
    
    # read in all SJ.tab files detected by 2-pass STAR and get the columns: chrom, intron start, intron end, strand, novelty
    junctionDFs = [pandas.read_table(inputFile, delim_whitespace = True, dtype = str, header = None, names = ['chrom', 'intron_start', 'intron_end', 'strandSTAR', 'motif', 'twoPassAnnotStatus'],
                                     usecols=[0, 1, 2, 3, 4, 5]) for inputFile in junctionFiles]
        
    # merge all SJ.tab files after column extraction and only keep unique ones
    junctions = pandas.concat(junctionDFs, ignore_index = True)
    junctions.drop_duplicates(keep = 'last', inplace = True)
    
    # gtf junctions
    known_juncs = pandas.read_table(annotation, delim_whitespace = True, dtype = str, header = None, names = ['chrom', 'intron_start', 'intron_end', 'strandSTAR', 'extraInfo'])
    known_juncs['annotStatus'] = 'annotated'
    '''
    rename columns values to be more meaningful:
        strand -> 0 = ud (undetermined), 1 = + strand, 2 = - strand
        motif -> 0 = nc (non-canonical), 1 = GT/AG, 2 = CT/AC, 3 = GC/AG, 4 = CT/GC, 5 = AT/AC, 6 = GT/AT
        twoPassAnnotStatus -> 0 = unannotated/novel, 1 = annotated based on 2nd pass.  Does not mean junction is novel or annotated relative to gtf
    '''
    
    junctions.replace({'strandSTAR': {'0':'ud' , '1': '+', '2':'-'}, 
                       'motif':{'0': 'nc', '1':'GT/AG', '2':'CT/AC', '3':'GC/AG', '4':'CT/GC', '5':'AT/AC', '6':'GT/AT'},
                       'twoPassAnnotStatus':{'0':'novel', '1':'annotated'}}, 
                       inplace = True)
    
    # merge in known gtf annotations into junctions detected by STAR across all samples
    fully_annotated = junctions.merge(known_juncs, how = 'left', on = ['chrom', 'intron_start', 'intron_end', 'strandSTAR'])
    fully_annotated['annotStatus'].fillna('novel', inplace=True)
    fully_annotated.drop(['twoPassAnnotStatus', 'extraInfo'], axis=1, inplace = True)
    
    return fully_annotated


@profile
def intron_retention(junctions : pandas.DataFrame, spladderOut : str, outdir : str) -> None:
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
    
    events_of_interest = filtering_criteria(event_type = 'intron_retention', base_df = info_combine, geneMatch = args.geneMatchFile, annot_filter = args.novel, diff_exp = args.deTestResults, pval_adj = args.pvalAdj)
    
    # free up memory
    del info_combine
    gc.collect()
    
    # only get unique event ids.-- drop all duplicated event ids and store counts of duplication in df for later annotation as warning for false positive of novel junction
    unique_events_of_interest = events_of_interest.drop(['annotStatus'], axis=1)
    unique_events_of_interest.drop_duplicates(keep = 'first', inplace = True)
    # total_events = unique_events_of_interest['event_id'].value_counts()
    '''
    Experimental blast implementation
    *** needs testing ***
    '''
    if args.frameMatch:
        
        unique_events_of_interest.reset_index(drop=True, inplace = True)

        unique_events_of_interest[['upstream_region_translated_frame', 'annotated_peptide', 'blast_expected_value', 'blast_gene']] = ''
        for idx, row in unique_events_of_interest.iterrows():
            print(idx)
            if row['strand'] == '+':
                blast_results = blast_search(blast = args.blast, db = args.blastDb, fasta = dna_fasta, chrom = row['chrom'], const_region_start = int(row['exon1_start']), const_region_end = int(row['exon1_end']), strand = row['strand'], gene = row['symbol'])
                #blast_results = blast_search(blast = blast, db = db, fasta = dna_fasta, chrom = row['chrom'], const_region_start = int(row['exon1_start']), const_region_end = int(row['exon1_end']), strand = row['strand'], gene = row['symbol'])

            elif row['strand'] == '-':
                blast_results = blast_search(blast = args.blast, db = args.blastDb, fasta = dna_fasta, chrom = row['chrom'], const_region_start = int(row['exon2_start']), const_region_end = int(row['exon2_end']), strand = row['strand'], gene = row['symbol'])
                #blast_results = blast_search(blast = blast, db = db, fasta = dna_fasta, chrom = row['chrom'], const_region_start = int(row['exon2_start']), const_region_end = int(row['exon2_end']), strand = row['strand'], gene = row['symbol'])

            #row[['upstream_region_translated_frame', 'annotated_peptide', 'blast_expected_value', 'blast_gene']] = blast_results[0], blast_results[1], blast_results[2], blast_results[3]
            unique_events_of_interest.at[unique_events_of_interest.index[idx], 'upstream_region_translated_frame'] = blast_results[0]
            unique_events_of_interest.at[unique_events_of_interest.index[idx], 'annotated_peptide'] = blast_results[1]
            unique_events_of_interest.at[unique_events_of_interest.index[idx], 'blast_expected_value'] = blast_results[2]
            unique_events_of_interest.at[unique_events_of_interest.index[idx], 'blast_gene'] = blast_results[3]

            del blast_results
            gc.collect()
    
    unique_events_of_interest.to_csv("testing_blast.txt", sep = "\t", index=False)
    filter_annotations = unique_events_of_interest.loc[unique_events_of_interest['annotated_peptide'] != 'none found'] # remove rows where a blast hit to a protein was not found
    unique_events_of_interest = filter_annotations
    del filter_annotations
    gc.collect()
    '''
    end of experimental blast implementation
    '''    
    
    # check if strandSTAR and strand are concordant
    #       if not concordant, check if strandSTAR is ud
    #            if ud, then use strand
    #       if concordant then continue   
    unique_events_of_interest.reset_index(drop=True, inplace = True)
    unique_events_of_interest['translated_ase_peptide'] = ''
    for idx,row in unique_events_of_interest.iterrows():
        if ((row['strandSTAR'] != row['strand']) & (row['strandSTAR'] == 'ud')):
            orfs, check_frame_region = calculate_orfs(event_type = args.eventType, fasta = dna_fasta, kmer_length = args.kmer, strand = row['strand'], chrom = row['chrom'], flank_left_start = int(row['exon1_start']), flank_left_end = int(row['exon1_end']), flank_right_start = int(row['exon2_start']), flank_right_end = int(row['exon2_end']), ase_start = int(row['intron_start']), ase_end = int(row['intron_end']))
            peptides = translate_orfs(event_id = row['event_id'], orfs =  orfs, peptide_bank = peptides, upstream_const = row['upstream_region_translated_frame'], correct_frame_only = args.frameMatch, frame_check = check_frame_region)
            unique_events_of_interest.at[unique_events_of_interest.index[idx], 'translated_ase_peptide'] = peptides[row['event_id']]
            print(peptides[row['event_id']])

        elif ((row['strandSTAR'] != row['strand']) & (row['strandSTAR'] != 'ud')):
            try:
                if args.testMode: # forces the translation despite mismatch or NaN and use spladder's strand since just for testing mode
                    orfs, check_frame_region = calculate_orfs(event_type = args.eventType, fasta = dna_fasta, kmer_length = args.kmer, strand = row['strand'], chrom = row['chrom'], flank_left_start = int(row['exon1_start']), flank_left_end = int(row['exon1_end']), flank_right_start = int(row['exon2_start']), flank_right_end = int(row['exon2_end']), ase_start = int(row['intron_start']), ase_end = int(row['intron_end']))
                    peptides = translate_orfs(event_id = row['event_id'], orfs = orfs, peptide_bank = peptides, upstream_const = row['upstream_region_translated_frame'], correct_frame_only = args.frameMatch, frame_check = check_frame_region)
                    unique_events_of_interest.at[unique_events_of_interest.index[idx], 'translated_ase_peptide'] = peptides[row['event_id']]
                    print(peptides[row['event_id']])

            except:
                print('star strand is {} and spladder strands is {}.  Skipping event.'.format(str(row['strandSTAR']), row['strand']))
                continue
            
        elif (row['strandSTAR'] == row['strand']):
            orfs, check_frame_region = calculate_orfs(event_type = args.eventType, fasta = dna_fasta, kmer_length = args.kmer, strand = row['strandSTAR'], chrom = row['chrom'], flank_left_start = int(row['exon1_start']), flank_left_end = int(row['exon1_end']), flank_right_start = int(row['exon2_start']), flank_right_end = int(row['exon2_end']), ase_start = int(row['intron_start']), ase_end = int(row['intron_end']))
            peptides = translate_orfs(event_id = row['event_id'], orfs = orfs, peptide_bank = peptides, upstream_const = row['upstream_region_translated_frame'], correct_frame_only = args.frameMatch, frame_check = check_frame_region)
            unique_events_of_interest.at[unique_events_of_interest.index[idx], 'translated_ase_peptide'] = peptides[row['event_id']]
            print(peptides[row['event_id']])

        del orfs
        gc.collect()

    # TO DO: calcuate percent overlap of kmer with flanking region
    
    
    # get psi info columns
    unique_events_of_interest.filter(regex = '.psi', axis = 1)
    # get read coverage info columns
    unique_events_of_interest.filter(regex = '_cov', axis = 1)
    # get mean event count and mean gene expression columns
    unique_events_of_interest.filter(regex = '^mean_', axis = 1)
    # get log2FC of differential testing in event counts and gene expression levels columns
    unique_events_of_interest.filter(regex = '^log2FC_', axis = 1)
   

def exon_skip(junctions : pandas.DataFrame, spladderOut : str, outdir : str) -> None:
    '''
    input:
    description:
    return:
    '''
    

    # dictionary that stores all ORF peptides for each event
    peptides = {}
    match_frame_results = [] # stores list of collect_results output from parallelizatoin of match_frame, assuming --frameMatch is set

    
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
    events_of_interest = filtering_criteria(event_type = 'exon_skip', base_df = final, geneMatch = args.geneMatchFile, annot_filter = args.novel, diff_exp = args.deTestResults, pval_adj = args.pvalAdj)

    # free up memory
    del tmp, final, info_combine
    gc.collect()
    
    # only get unique event ids.-- drop all duplicated event ids and store counts of duplication in df for later annotation as warning for false positive of novel junction
    unique_events_of_interest = events_of_interest.drop(['annotStatus_pre', 'annotStatus_aft'], axis=1)
    unique_events_of_interest.drop_duplicates(keep = 'first', inplace = True)
    total_events = unique_events_of_interest['event_id'].value_counts()

    # function to collect results of match_frame, when --frameMatch flag is set; used for parallelization coordination
    #def collect_results(result:pandas.DataFrame) -> None:
    #    match_frame_results.extend(result)
    
    def match_frame(chunk_df : pandas.DataFrame, blast_tmp_file : int, fasta : Bio.File._IndexedSeqFileDict, shared_obj : list) -> list:
        '''
        Experimental blast implementation
        *** needs testing ***
        '''
        # each process requires its own instance of fasta
        fasta = SeqIO.index(args.fasta, 'fasta')
        
        chunk_df.reset_index(drop=True, inplace = True)

        chunk_df[['upstream_region_translated_frame', 'annotated_peptide', 'blast_expected_value', 'blast_gene']] = ''
        for idx, row in chunk_df.iterrows():
            print(idx)
            print(row['event_id'], row['chrom'], row['exon_pre_start'], row['exon_pre_end'], row['symbol'])
            if row['strand'] == '+':
                blast_results = blast_search(event_id = row['event_id'], blast = args.blast, db = args.blastDb, fasta = fasta, chrom = row['chrom'], const_region_start = int(row['exon_pre_start']), const_region_end = int(row['exon_pre_end']), strand = row['strand'], gene = row['symbol'], blast_tmp_file = blast_tmp_file)
                #blast_results = blast_search(blast = blast, db = db, fasta = dna_fasta, chrom = row['chrom'], const_region_start = int(row['exon_pre_start']), const_region_end = int(row['exon_pre_end']), strand = row['strand'], gene = row['symbol_x'])

            elif row['strand'] == '-':
                blast_results = blast_search(event_id = row['event_id'], blast = args.blast, db = args.blastDb, fasta = fasta, chrom = row['chrom'], const_region_start = int(row['exon_aft_start']), const_region_end = int(row['exon_aft_end']), strand = row['strand'], gene = row['symbol'], blast_tmp_file = blast_tmp_file)
                #blast_results = blast_search(blast = blast, db = db, fasta = dna_fasta, chrom = row['chrom'], const_region_start = int(row['exon_aft_start']), const_region_end = int(row['exon_aft_end']), strand = row['strand'], gene = row['symbol_x'])

            row[['upstream_region_translated_frame', 'annotated_peptide', 'blast_expected_value', 'blast_gene']] = blast_results[0], blast_results[1], blast_results[2], blast_results[3]
            chunk_df.at[chunk_df.index[idx], 'upstream_region_translated_frame'] = blast_results[0]
            chunk_df.at[chunk_df.index[idx], 'annotated_peptide'] = blast_results[1]
            chunk_df.at[chunk_df.index[idx], 'blast_expected_value'] = blast_results[2]
            chunk_df.at[chunk_df.index[idx], 'blast_gene'] = blast_results[3]

            del blast_results
            gc.collect()
        
        shared_obj.append(chunk_df)
        return shared_obj # converts pandas DF to list
            
    if args.frameMatch:
        
        number_of_processes = math.ceil(multiprocessing.cpu_count()*.70) # use 70% of available CPUs

        df_chunks = numpy.array_split(unique_events_of_interest, number_of_processes) # returns a list of dataframes split into chunks matching number of processes available
        
        shared_list = multiprocessing.Manager().list()
        
        parallelized_jobs = [multiprocessing.Process(target=match_frame, args = (df_chunks[dfs], dfs, dna_fasta, shared_list, )) for dfs in range(0, len(df_chunks))]
        
        for processes in parallelized_jobs:
            processes.start()
        
        for processes in parallelized_jobs:
            processes.join()
        
        unique_events_of_interest = pandas.concat(shared_list)
        
        unique_events_of_interest.to_csv("test_parallel_exon_skip_testing_blast_03212021.txt", sep = '\t', index=False)
    
    '''
    TODO add --outdir args for any read and write csv functions
    '''
    unique_events_of_interest = pandas.read_csv("test_parallel_exon_skip_testing_blast_03212021.txt", sep = '\t')
    filter_annotations = unique_events_of_interest.loc[unique_events_of_interest['annotated_peptide'] != 'none found'] # remove rows where a blast hit to a protein was not found
    unique_events_of_interest = filter_annotations
    del filter_annotations
    gc.collect()

    # check if strandSTAR and strand are concordant
    #       if not concordant, check if strandSTAR is ud
    #            if ud, then use strand
    #       if concordant then continue
    unique_events_of_interest.reset_index(drop=True, inplace = True)
    unique_events_of_interest['translated_ase_peptide'] = ''
        
    ### parallelization of calculation and tranlation of orfs
    
    def get_peptides(chunk_df : pandas.DataFrame, shared_peptides : Dict[str, Dict[str,str]], shared_dfs : list ) -> Tuple[Dict[str, Dict[str,str]], list]:    
        
        fasta = SeqIO.index(args.fasta, 'fasta')
        
        # reset index so that each dataframe chunk has the right index ID corresponding to the iterrows() index
        chunk_df.reset_index(drop=True, inplace = True)
        
        for idx,row in chunk_df.iterrows():
            print(idx)
            try:
                assert row['strandSTAR_pre'] == row['strandSTAR_aft'], 'There is a strand conflict at the following event_id: {}.  Skipping this event.'.format(row['event_id'])
            except  AssertionError as err:
                print(err)
                continue
            
            # args.kmer has one subtracted since the ase event is built into the upsteam and downstream contant exons    
            
            if ((row['strandSTAR_pre'] != row['strand']) & (row['strandSTAR_pre'] == 'ud')): # note it does not matter if we use strandSTAR_pre for strandStar_aft since we already checked the assertion that the STAR strands are concordant with each other
                #orfs = calculate_orfs(event_type = args.eventType, fasta = dna_fasta, kmer_length = args.kmer, strand = row['strand'], chrom = row['chrom'], flank_left_start = int(row['exon_pre_start']), flank_left_end = int(row['exon_pre_end']), flank_right_start = int(row['exon_aft_start']), flank_right_end = int(row['exon_aft_end']), ase_start = int(row['exon_start']), ase_end = int(row['exon_end']))
                orfs, check_frame_region = calculate_orfs(event_type = args.eventType, fasta = fasta, kmer_length = args.kmer - 1, strand = row['strand'], chrom = row['chrom'], flank_left_start = int(row['exon_pre_start']), flank_left_end = int(row['exon_pre_end'])-3, flank_right_start = int(row['exon_aft_start'])+3, flank_right_end = int(row['exon_aft_end']), ase_start = int(row['exon_pre_end'])-3, ase_end = int(row['exon_aft_start'])+3)
                shared_peptides = translate_orfs(event_id = row['event_id'], orfs =  orfs, peptide_bank = shared_peptides, upstream_const = row['upstream_region_translated_frame'], correct_frame_only = args.frameMatch, frame_check = check_frame_region)
                chunk_df.at[chunk_df.index[idx], 'translated_ase_peptide'] = shared_peptides[row['event_id']]

            elif ((row['strandSTAR_pre'] != row['strand']) & (row['strandSTAR_pre'] != 'ud')):
                try:
                    if args.testMode: # forces the translation despite mismatch or NaN and use spladder's strand since just for testing mode
                        orfs, check_frame_region = calculate_orfs(event_type = args.eventType, fasta = fasta, kmer_length = args.kmer - 1, strand = row['strand'], chrom = row['chrom'], flank_left_start = int(row['exon_pre_start']), flank_left_end = int(row['exon_pre_end'])-3, flank_right_start = int(row['exon_aft_start'])+3, flank_right_end = int(row['exon_aft_end']), ase_start = int(row['exon_pre_end'])-3, ase_end = int(row['exon_aft_start'])+3)
                        shared_peptides = translate_orfs(event_id = row['event_id'], orfs = orfs, peptide_bank = shared_peptides, upstream_const = row['upstream_region_translated_frame'], correct_frame_only = args.frameMatch, frame_check = check_frame_region)
                        chunk_df.at[chunk_df.index[idx], 'translated_ase_peptide'] = shared_peptides[row['event_id']]

                except:
                        print('star strand is {} and spladder strands is {}.  Skipping event.'.format(row['strandSTAR_pre'], row['strand']))
                        continue
            
            elif (row['strandSTAR_pre'] == row['strand']):
                orfs, check_frame_region = calculate_orfs(event_type = args.eventType, fasta = fasta, kmer_length = args.kmer - 1, strand = row['strandSTAR_pre'], chrom = row['chrom'], flank_left_start = int(row['exon_pre_start']), flank_left_end = int(row['exon_pre_end'])-3, flank_right_start = int(row['exon_aft_start'])+3, flank_right_end = int(row['exon_aft_end']), ase_start = int(row['exon_pre_end'])-3, ase_end = int(row['exon_aft_start'])+3)
                shared_peptides = translate_orfs(event_id = row['event_id'], orfs = orfs, peptide_bank = shared_peptides, upstream_const = row['upstream_region_translated_frame'], correct_frame_only = args.frameMatch, frame_check = check_frame_region)
                chunk_df.at[chunk_df.index[idx], 'translated_ase_peptide'] = shared_peptides[row['event_id']]

            del orfs
            gc.collect() 
        
        shared_dfs.append(chunk_df)

        return shared_peptides, shared_dfs
        
    number_of_processes = math.ceil(multiprocessing.cpu_count()*.70) # use 70% of available CPUs
    df_chunks = numpy.array_split(unique_events_of_interest, number_of_processes) # returns a list of dataframes split into chunks matching number of processes available
    shared_peptides = multiprocessing.Manager().dict()
    shared_dfs = multiprocessing.Manager().list()  
    parallelized_jobs = [multiprocessing.Process(target=get_peptides, args = (df_chunks[dfs], shared_peptides, shared_dfs, )) for dfs in range(0, len(df_chunks))]
        
    for processes in parallelized_jobs:
        processes.start()
        
    for processes in parallelized_jobs:
        processes.join()
    
    unique_events_of_interest = pandas.concat(shared_dfs)
        
    unique_events_of_interest.to_csv("first_10_test_parallel_exon_skip_testing_blast_03292022_RESULTS.txt", sep = '\t', index=False)


    '''

    # TO DO: calcuate percent overlap of kmer with flanking region
    
    unique_events_of_interest.to_csv("exon_skip_final_dataframe_test.txt", sep = "\t")
    # get psi info columns
    unique_events_of_interest.filter(regex = '.psi', axis = 1)
    # get read coverage info columns
    unique_events_of_interest.filter(regex = '_cov', axis = 1)
    # get mean event count and mean gene expression columns
    unique_events_of_interest.filter(regex = '^mean_', axis = 1)
    # get log2FC of differential testing in event counts and gene expression levels columns
    unique_events_of_interest.filter(regex = '^log2FC_', axis = 1)

   '''


def three_prime_alt(junctions : pandas.DataFrame, spladderOut : str, outdir : str) -> None:
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

 
    rev_strand_events['alt1_intron_start'] = rev_strand_events['exon_alt1_end'].astype(int) + 1
    rev_strand_events['alt1_intron_end'] = rev_strand_events['exon_const_start'].astype(int) - 1
    rev_strand_events['alt2_intron_start'] = rev_strand_events['exon_alt2_end'].astype(int) + 1
    rev_strand_events['alt2_intron_end'] = rev_strand_events['exon_const_start'].astype(int) - 1
    
    
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
    fwd_alt1_info_combine = fwd_strand_events.merge(fwd_tmp, how = 'inner', on = ['chrom', 'alt1_intron_start', 'alt1_intron_end'])
    
    fwd_tmp = junctions.rename({'intron_start':'alt2_intron_start', 'intron_end':'alt2_intron_end', 'strandSTAR':'strandSTAR_alt2', 'motif':'motif_alt2', 'annotStatus':'annotStatus_alt2' }, axis = 'columns')
    fwd_final = fwd_alt1_info_combine.merge(fwd_tmp, how = 'left', on = ['chrom', 'alt2_intron_start', 'alt2_intron_end'])

      
    # merge star intron annotation into reverse strand events
    rev_tmp = junctions.rename({'intron_start':'alt1_intron_start', 'intron_end':'alt1_intron_end', 'strandSTAR':'strandSTAR_alt1', 'motif':'motif_alt1', 'annotStatus':'annotStatus_alt1' }, axis = 'columns')
    rev_alt1_info_combine = rev_strand_events.merge(rev_tmp, how = 'left', on = ['chrom', 'alt1_intron_start', 'alt1_intron_end'])

    rev_tmp = junctions.rename({'intron_start':'alt2_intron_start', 'intron_end':'alt2_intron_end', 'strandSTAR':'strandSTAR_alt2', 'motif':'motif_alt2', 'annotStatus':'annotStatus_alt2' }, axis = 'columns')
    rev_final = rev_alt1_info_combine.merge(rev_tmp, how = 'left', on = ['chrom', 'alt2_intron_start', 'alt2_intron_end'])
    

    # merge dataframes back together
    combined_strand_events = pandas.concat([fwd_final, rev_final], axis = 0, ignore_index= True)
    
    # free memory
    del fwd_tmp, fwd_alt1_info_combine, rev_tmp, rev_alt1_info_combine, fwd_final, rev_final
    gc.collect()
    
    # apply data filtering criteria
    events_of_interest = filtering_criteria(event_type = 'three_prime_alt', base_df = combined_strand_events, geneMatch = args.geneMatchFile, annot_filter = args.novel, diff_exp = args.deTestResults, pval_adj = args.pvalAdj)

    # only get unique event ids.-- drop all duplicated event ids and store counts of duplication in df for later annotation as warning for false positive of novel junction
    unique_events_of_interest = events_of_interest.drop(['annotStatus_alt1', 'annotStatus_alt2'], axis=1)
    unique_events_of_interest.drop_duplicates(keep = 'first', inplace = True)
    total_events = unique_events_of_interest['event_id'].value_counts()

    '''
    Experimental blast implementation
    *** needs testing ***
    '''
    
    # TO DO
    
    '''
    end of experimental blast implementation
    '''  
    
    # check if strandSTAR and strand are concordant
    #       if alt1 and alt2 are not concordant within just star, then skip\
    #       if not concordant between star and stran, check if strandSTAR is ud
    #            if ud, then use strand
    #       if concordant then continue

    for idx,row in unique_events_of_interest.iterrows():
        try:
            assert row['strandSTAR_alt1'] == row['strandSTAR_alt2'], 'There is a strand conflict at the following event_id: {}.  Skipping this event.'.format(row['event_id'])
        except  AssertionError as err:
            print(err)
            continue
            
        # check strand for alt events since genomic positions are in reverse order for this alternative splicing event type
        
        if ((row['strandSTAR_alt1'] != row['strand']) & (row['strandSTAR_alt1'] == 'ud')): # note it does not matter if we use strandSTAR_pre for strandStar_aft since we already checked the assertion that the STAR strands are concordant with each other
            if row['strand'] == '+':
                orfs_alt1 = calculate_orfs(event_type = args.eventType, fasta = dna_fasta, kmer_length = args.kmer, strand = row['strand'], chrom = row['chrom'], flank_left_start = int(row['exon_const_start']), flank_left_end = int(row['exon_const_end']), flank_right_start = int(row['exon_alt1_start']), flank_right_end = int(row['exon_alt1_end']), ase_start = int(row['exon_alt1_start']), ase_end = None)
                peptides_alt1 = translate_orfs(event_id = row['event_id'], orfs =  orfs_alt1, peptide_bank = peptides_alt1)
                del orfs_alt1
                gc.collect()
                orfs_alt2 = calculate_orfs(event_type = args.eventType, fasta = dna_fasta, kmer_length = args.kmer, strand = row['strand'], chrom = row['chrom'], flank_left_start = int(row['exon_const_start']), flank_left_end = int(row['exon_const_end']), flank_right_start = int(row['exon_alt2_start']), flank_right_end = int(row['exon_alt2_end']), ase_start = int(row['exon_alt2_start']), ase_end = None)
                peptides_alt2 = translate_orfs(event_id = row['event_id'], orfs =  orfs_alt2, peptide_bank = peptides_alt2)
                del orfs_alt2
                gc.collect()
                
                
            elif row['strand'] =='-':
                orfs_alt1 = calculate_orfs(event_type = args.eventType, fasta = dna_fasta, kmer_length = args.kmer, strand = row['strand'], chrom = row['chrom'], flank_left_start = int(row['exon_alt1_start']), flank_left_end = int(row['exon_alt1_end']), flank_right_start = int(row['exon_const_start']), flank_right_end = int(row['exon_const_end']), ase_start = int(row['exon_alt1_end']), ase_end = None)
                peptides_alt1 = translate_orfs(event_id = row['event_id'], orfs =  orfs_alt1, peptide_bank = peptides_alt1)
                del orfs_alt1
                gc.collect()
                orfs_alt2 = calculate_orfs(event_type = args.eventType, fasta = dna_fasta, kmer_length = args.kmer, strand = row['strand'], chrom = row['chrom'], flank_left_start = int(row['exon_alt2_start']), flank_left_end = int(row['exon_alt2_end']), flank_right_start = int(row['exon_const_start']), flank_right_end = int(row['exon_const_end']), ase_start = int(row['exon_alt2_end']), ase_end = None)
                peptides_alt2 = translate_orfs(event_id = row['event_id'], orfs =  orfs_alt2, peptide_bank = peptides_alt2)
                del orfs_alt2
                gc.collect()


        elif ((row['strandSTAR_alt1'] != row['strand']) & (row['strandSTAR_alt1'] != 'ud')):
            try:
                if args.testMode: # forces the translation despite mismatch or NaN and use spladder's strand since just for testing mode
                    if row['strand'] == '+':
                        orfs_alt1 = calculate_orfs(event_type = args.eventType, fasta = dna_fasta, kmer_length = args.kmer, strand = row['strand'], chrom = row['chrom'], flank_left_start = int(row['exon_const_start']), flank_left_end = int(row['exon_const_end']), flank_right_start = int(row['exon_alt1_start']), flank_right_end = int(row['exon_alt1_end']), ase_start = int(row['exon_alt1_start']), ase_end = None)
                        peptides_alt1 = translate_orfs(event_id = row['event_id'], orfs =  orfs_alt1, peptide_bank = peptides_alt1)
                        del orfs_alt1
                        gc.collect()
                        orfs_alt2 = calculate_orfs(event_type = args.eventType, fasta = dna_fasta, kmer_length = args.kmer, strand = row['strand'], chrom = row['chrom'], flank_left_start = int(row['exon_const_start']), flank_left_end = int(row['exon_const_end']), flank_right_start = int(row['exon_alt2_start']), flank_right_end = int(row['exon_alt2_end']), ase_start = int(row['exon_alt2_start']), ase_end = None)
                        peptides_alt2 = translate_orfs(event_id = row['event_id'], orfs =  orfs_alt2, peptide_bank = peptides_alt2)
                        del orfs_alt2
                        gc.collect()

                    elif row['strand'] =='-':
                        orfs_alt1 = calculate_orfs(event_type = args.eventType, fasta = dna_fasta, kmer_length = args.kmer, strand = row['strand'], chrom = row['chrom'], flank_left_start = int(row['exon_alt1_start']), flank_left_end = int(row['exon_alt1_end']), flank_right_start = int(row['exon_const_start']), flank_right_end = int(row['exon_const_end']), ase_start = int(row['exon_alt1_end']), ase_end = None)
                        peptides_alt1 = translate_orfs(event_id = row['event_id'], orfs =  orfs_alt1, peptide_bank = peptides_alt1)
                        del orfs_alt1
                        gc.collect()
                        orfs_alt2 = calculate_orfs(event_type = args.eventType, fasta = dna_fasta, kmer_length = args.kmer, strand = row['strand'], chrom = row['chrom'], flank_left_start = int(row['exon_alt2_start']), flank_left_end = int(row['exon_alt2_end']), flank_right_start = int(row['exon_const_start']), flank_right_end = int(row['exon_const_end']), ase_start = int(row['exon_alt2_end']), ase_end = None)
                        peptides_alt2 = translate_orfs(event_id = row['event_id'], orfs =  orfs_alt2, peptide_bank = peptides_alt2)
                        del orfs_alt2
                        gc.collect()
            except:
                    print('star strand is {} and spladder strands is {}.  Skipping event.'.format(row['strandSTAR_pre'], row['strand']))
                    continue
        
        
        elif (row['strandSTAR_alt1'] == row['strand']):
            
            if row['strand'] == '+':
                orfs_alt1 = calculate_orfs(event_type = args.eventType, fasta = dna_fasta, kmer_length = args.kmer, strand = row['strandSTAR_alt1'], chrom = row['chrom'], flank_left_start = int(row['exon_const_start']), flank_left_end = int(row['exon_const_end']), flank_right_start = int(row['exon_alt1_start']), flank_right_end = int(row['exon_alt1_end']), ase_start = int(row['exon_alt1_start']), ase_end = None)
                peptides_alt1 = translate_orfs(event_id = row['event_id'], orfs =  orfs_alt1, peptide_bank = peptides_alt1)
                del orfs_alt1
                gc.collect()
                orfs_alt2 = calculate_orfs(event_type = args.eventType, fasta = dna_fasta, kmer_length = args.kmer, strand = row['strandSTAR_alt2'], chrom = row['chrom'], flank_left_start = int(row['exon_const_start']), flank_left_end = int(row['exon_const_end']), flank_right_start = int(row['exon_alt2_start']), flank_right_end = int(row['exon_alt2_end']), ase_start = int(row['exon_alt2_start']), ase_end = None)
                peptides_alt2 = translate_orfs(event_id = row['event_id'], orfs =  orfs_alt2, peptide_bank = peptides_alt2)
                del orfs_alt2
                gc.collect()

            elif row['strand'] =='-':
                orfs_alt1 = calculate_orfs(event_type = args.eventType, fasta = dna_fasta, kmer_length = args.kmer, strand = row['strandSTAR_alt1'], chrom = row['chrom'], flank_left_start = int(row['exon_alt1_start']), flank_left_end = int(row['exon_alt1_end']), flank_right_start = int(row['exon_const_start']), flank_right_end = int(row['exon_const_end']), ase_start = int(row['exon_alt1_end']), ase_end = None)
                peptides_alt1 = translate_orfs(event_id = row['event_id'], orfs =  orfs_alt1, peptide_bank = peptides_alt1)
                del orfs_alt1
                gc.collect()
                orfs_alt2 = calculate_orfs(event_type = args.eventType, fasta = dna_fasta, kmer_length = args.kmer, strand = row['strandSTAR_alt2'], chrom = row['chrom'], flank_left_start = int(row['exon_alt2_start']), flank_left_end = int(row['exon_alt2_end']), flank_right_start = int(row['exon_const_start']), flank_right_end = int(row['exon_const_end']), ase_start = int(row['exon_alt2_end']), ase_end = None)
                peptides_alt2 = translate_orfs(event_id = row['event_id'], orfs =  orfs_alt2, peptide_bank = peptides_alt2)
                del orfs_alt2
                gc.collect()
   
    # TO DO: calcuate percent overlap of kmer with flanking region
    
    
    # get psi info columns
    unique_events_of_interest.filter(regex = '.psi', axis = 1)
    # get read coverage info columns
    unique_events_of_interest.filter(regex = '_cov', axis = 1)
    # get mean event count and mean gene expression columns
    unique_events_of_interest.filter(regex = '^mean_', axis = 1)
    # get log2FC of differential testing in event counts and gene expression levels columns
    unique_events_of_interest.filter(regex = '^log2FC_', axis = 1)
       
    
    
def five_prime_alt(junctions : pandas.DataFrame, spladderOut : str, outdir : str) -> None:
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
    

    fwd_strand_events['alt1_intron_start'] = fwd_strand_events['exon_alt1_end'].astype(int) + 1
    fwd_strand_events['alt1_intron_end'] = fwd_strand_events['exon_const_start'].astype(int) - 1
    fwd_strand_events['alt2_intron_start'] = fwd_strand_events['exon_alt2_end'].astype(int) + 1
    fwd_strand_events['alt2_intron_end'] = fwd_strand_events['exon_const_start'].astype(int) - 1

 
    rev_strand_events['alt1_intron_start'] = rev_strand_events['exon_const_end'].astype(int) + 1
    rev_strand_events['alt1_intron_end'] = rev_strand_events['exon_alt1_start'].astype(int) - 1
    rev_strand_events['alt2_intron_start'] = rev_strand_events['exon_const_end'].astype(int) + 1
    rev_strand_events['alt2_intron_end'] = rev_strand_events['exon_alt2_start'].astype(int) - 1
    
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
    fwd_alt1_info_combine = fwd_strand_events.merge(fwd_tmp, how = 'inner', on = ['chrom', 'alt1_intron_start', 'alt1_intron_end'])
    
    fwd_tmp = junctions.rename({'intron_start':'alt2_intron_start', 'intron_end':'alt2_intron_end', 'strandSTAR':'strandSTAR_alt2', 'motif':'motif_alt2', 'annotStatus':'annotStatus_alt2' }, axis = 'columns')
    fwd_final = fwd_alt1_info_combine.merge(fwd_tmp, how = 'left', on = ['chrom', 'alt2_intron_start', 'alt2_intron_end'])

      
    # merge star intron annotation into reverse strand events
    rev_tmp = junctions.rename({'intron_start':'alt1_intron_start', 'intron_end':'alt1_intron_end', 'strandSTAR':'strandSTAR_alt1', 'motif':'motif_alt1', 'annotStatus':'annotStatus_alt1' }, axis = 'columns')
    rev_alt1_info_combine = rev_strand_events.merge(rev_tmp, how = 'left', on = ['chrom', 'alt1_intron_start', 'alt1_intron_end'])

    rev_tmp = junctions.rename({'intron_start':'alt2_intron_start', 'intron_end':'alt2_intron_end', 'strandSTAR':'strandSTAR_alt2', 'motif':'motif_alt2', 'annotStatus':'annotStatus_alt2' }, axis = 'columns')
    rev_final = rev_alt1_info_combine.merge(rev_tmp, how = 'left', on = ['chrom', 'alt2_intron_start', 'alt2_intron_end'])
    
    # merge dataframes back together
    combined_strand_events = pandas.concat([fwd_final, rev_final], axis = 0, ignore_index= True)
    
    # free memory
    del fwd_tmp, fwd_alt1_info_combine, rev_tmp, rev_alt1_info_combine, fwd_final, rev_final
    gc.collect()
    
    # apply data filtering criteria
    events_of_interest = filtering_criteria(event_type = 'five_prime_alt', base_df = combined_strand_events, geneMatch = args.geneMatchFile, annot_filter = args.novel, diff_exp = args.deTestResults, pval_adj = args.pvalAdj)

    # only get unique event ids.-- drop all duplicated event ids and store counts of duplication in df for later annotation as warning for false positive of novel junction
    unique_events_of_interest = events_of_interest.drop(['annotStatus_alt1', 'annotStatus_alt2'], axis=1)
    unique_events_of_interest.drop_duplicates(keep = 'first', inplace = True)
    total_events = unique_events_of_interest['event_id'].value_counts()
    
    
    '''
    Experimental blast implementation
    *** needs testing ***
    '''
    
    # TO DO
    
    '''
    end of experimental blast implementation
    '''  
   
   
    # check if strandSTAR and strand are concordant
    #       if alt1 and alt2 are not concordant within just star, then skip\
    #       if not concordant between star and strand, check if strandSTAR is ud
    #            if ud, then use strand
    #       if concordant then continue

    for idx,row in unique_events_of_interest.iterrows():
        try:
            assert row['strandSTAR_alt1'] == row['strandSTAR_alt2'], 'There is a strand conflict at the following event_id: {}.  Skipping this event.'.format(row['event_id'])
        except  AssertionError as err:
            print(err)
            continue
            
        # check strand for alt events since genomic positions are in reverse order for this alternative splicing event type
        
        if ((row['strandSTAR_alt1'] != row['strand']) & (row['strandSTAR_alt1'] == 'ud')): # note it does not matter if we use strandSTAR_pre for strandStar_aft since we already checked the assertion that the STAR strands are concordant with each other
            if row['strand'] == '+':
                orfs_alt1 = calculate_orfs(event_type = args.eventType, fasta = dna_fasta, kmer_length = args.kmer, strand = row['strand'], chrom = row['chrom'], flank_left_start = int(row['exon_alt1_start']), flank_left_end = int(row['exon_alt1_end']), flank_right_start = int(row['exon_const_start']), flank_right_end = int(row['exon_const_end']), ase_start = int(row['exon_alt1_end']), ase_end = None)
                peptides_alt1 = translate_orfs(event_id = row['event_id'], orfs =  orfs_alt1, peptide_bank = peptides_alt1)
                del orfs_alt1
                gc.collect()
                orfs_alt2 = calculate_orfs(event_type = args.eventType, fasta = dna_fasta, kmer_length = args.kmer, strand = row['strand'], chrom = row['chrom'], flank_left_start = int(row['exon_alt2_start']), flank_left_end = int(row['exon_alt2_end']), flank_right_start = int(row['exon_const_start']), flank_right_end = int(row['exon_const_end']), ase_start = int(row['exon_alt2_end']), ase_end = None)
                peptides_alt2 = translate_orfs(event_id = row['event_id'], orfs =  orfs_alt2, peptide_bank = peptides_alt2)
                del orfs_alt2
                gc.collect()                
            
            elif row['strand'] =='-':
                orfs_alt1 = calculate_orfs(event_type = args.eventType, fasta = dna_fasta, kmer_length = args.kmer, strand = row['strand'], chrom = row['chrom'], flank_left_start = int(row['exon_const_start']), flank_left_end = int(row['exon_const_end']), flank_right_start = int(row['exon_alt1_start']), flank_right_end = int(row['exon_alt1_end']), ase_start = int(row['exon_alt1_start']), ase_end = None)
                peptides_alt1 = translate_orfs(event_id = row['event_id'], orfs =  orfs_alt1, peptide_bank = peptides_alt1)
                del orfs_alt1
                gc.collect()
                orfs_alt2 = calculate_orfs(event_type = args.eventType, fasta = dna_fasta, kmer_length = args.kmer, strand = row['strand'], chrom = row['chrom'], flank_left_start = int(row['exon_const_start']), flank_left_end = int(row['exon_const_end']), flank_right_start = int(row['exon_alt2_start']), flank_right_end = int(row['exon_alt2_end']), ase_start = int(row['exon_alt2_start']), ase_end = None)
                peptides_alt2 = translate_orfs(event_id = row['event_id'], orfs =  orfs_alt2, peptide_bank = peptides_alt2)
                del orfs_alt2
                gc.collect()



        elif ((row['strandSTAR_alt1'] != row['strand']) & (row['strandSTAR_alt1'] != 'ud')):
            try:
                if args.testMode: # forces the translation despite mismatch or NaN and use spladder's strand since just for testing mode
                    if row['strand'] == '+':
                        orfs_alt1 = calculate_orfs(event_type = args.eventType, fasta = dna_fasta, kmer_length = args.kmer, strand = row['strand'], chrom = row['chrom'], flank_left_start = int(row['exon_alt1_start']), flank_left_end = int(row['exon_alt1_end']), flank_right_start = int(row['exon_const_start']), flank_right_end = int(row['exon_const_end']), ase_start = int(row['exon_alt1_end']), ase_end = None)
                        peptides_alt1 = translate_orfs(event_id = row['event_id'], orfs =  orfs_alt1, peptide_bank = peptides_alt1)
                        del orfs_alt1
                        gc.collect()
                        orfs_alt2 = calculate_orfs(event_type = args.eventType, fasta = dna_fasta, kmer_length = args.kmer, strand = row['strand'], chrom = row['chrom'], flank_left_start = int(row['exon_alt2_start']), flank_left_end = int(row['exon_alt2_end']), flank_right_start = int(row['exon_const_start']), flank_right_end = int(row['exon_const_end']), ase_start = int(row['exon_alt2_end']), ase_end = None)
                        peptides_alt2 = translate_orfs(event_id = row['event_id'], orfs =  orfs_alt2, peptide_bank = peptides_alt2)
                        del orfs_alt2
                        gc.collect()
    
                    elif row['strand'] =='-':
                        orfs_alt1 = calculate_orfs(event_type = args.eventType, fasta = dna_fasta, kmer_length = args.kmer, strand = row['strand'], chrom = row['chrom'], flank_left_start = int(row['exon_const_start']), flank_left_end = int(row['exon_const_end']), flank_right_start = int(row['exon_alt1_start']), flank_right_end = int(row['exon_alt1_end']), ase_start = int(row['exon_alt1_start']), ase_end = None)
                        peptides_alt1 = translate_orfs(event_id = row['event_id'], orfs =  orfs_alt1, peptide_bank = peptides_alt1)
                        del orfs_alt1
                        gc.collect()
                        orfs_alt2 = calculate_orfs(event_type = args.eventType, fasta = dna_fasta, kmer_length = args.kmer, strand = row['strand'], chrom = row['chrom'], flank_left_start = int(row['exon_const_start']), flank_left_end = int(row['exon_const_end']), flank_right_start = int(row['exon_alt2_start']), flank_right_end = int(row['exon_alt2_end']), ase_start = int(row['exon_alt2_start']), ase_end = None)
                        peptides_alt2 = translate_orfs(event_id = row['event_id'], orfs =  orfs_alt2, peptide_bank = peptides_alt2)
                        del orfs_alt2
                        gc.collect()
            except:
                print('star strand is {} and spladder strands is {}.  Skipping event.'.format(row['strandSTAR_alt1'], row['strand']))
                continue
        
        
        elif (row['strandSTAR_alt1'] == row['strand']):
            
            if row['strand'] == '+':
                orfs_alt1 = calculate_orfs(event_type = args.eventType, fasta = dna_fasta, kmer_length = args.kmer, strand = row['strandSTAR_alt1'], chrom = row['chrom'], flank_left_start = int(row['exon_alt1_start']), flank_left_end = int(row['exon_alt1_end']), flank_right_start = int(row['exon_const_start']), flank_right_end = int(row['exon_const_end']), ase_start = int(row['exon_alt1_end']), ase_end = None)
                peptides_alt1 = translate_orfs(event_id = row['event_id'], orfs =  orfs_alt1, peptide_bank = peptides_alt1)
                del orfs_alt1
                gc.collect()
                orfs_alt2 = calculate_orfs(event_type = args.eventType, fasta = dna_fasta, kmer_length = args.kmer, strand = row['strandSTAR_alt2'], chrom = row['chrom'], flank_left_start = int(row['exon_alt2_start']), flank_left_end = int(row['exon_alt2_end']), flank_right_start = int(row['exon_const_start']), flank_right_end = int(row['exon_const_end']), ase_start = int(row['exon_alt2_end']), ase_end = None)
                peptides_alt2 = translate_orfs(event_id = row['event_id'], orfs =  orfs_alt2, peptide_bank = peptides_alt2)
                del orfs_alt2
                gc.collect()
    
            elif row['strand'] =='-':
                orfs_alt1 = calculate_orfs(event_type = args.eventType, fasta = dna_fasta, kmer_length = args.kmer, strand = row['strandSTAR_alt1'], chrom = row['chrom'], flank_left_start = int(row['exon_const_start']), flank_left_end = int(row['exon_const_end']), flank_right_start = int(row['exon_alt1_start']), flank_right_end = int(row['exon_alt1_end']), ase_start = int(row['exon_alt1_start']), ase_end = None)
                peptides_alt1 = translate_orfs(event_id = row['event_id'], orfs =  orfs_alt1, peptide_bank = peptides_alt1)
                del orfs_alt1
                gc.collect()
                orfs_alt2 = calculate_orfs(event_type = args.eventType, fasta = dna_fasta, kmer_length = args.kmer, strand = row['strandSTAR_alt2'], chrom = row['chrom'], flank_left_start = int(row['exon_const_start']), flank_left_end = int(row['exon_const_end']), flank_right_start = int(row['exon_alt2_start']), flank_right_end = int(row['exon_alt2_end']), ase_start = int(row['exon_alt2_start']), ase_end = None)
                peptides_alt2 = translate_orfs(event_id = row['event_id'], orfs =  orfs_alt2, peptide_bank = peptides_alt2)
                del orfs_alt2
                gc.collect()
    
    
    # TO DO: calcuate percent overlap of kmer with flanking region
    
    # get psi info columns
    unique_events_of_interest.filter(regex = '.psi', axis = 1)
    # get read coverage info columns
    unique_events_of_interest.filter(regex = '_cov', axis = 1)
    # get mean event count and mean gene expression columns
    unique_events_of_interest.filter(regex = '^mean_', axis = 1)
    # get log2FC of differential testing in event counts and gene expression levels columns
    unique_events_of_interest.filter(regex = '^log2FC_', axis = 1)

    

def mutex(junctions : pandas.DataFrame, spladderOut : str, outdir : str) -> None:
    '''
    input:
    description:
    return:
    '''
    # dictionary that stores all ORF peptides for each event
    peptides_exon1 = {}
    peptides_exon2 = {}
    
    # read in data from spladder and then rename columns for df merge
    as_events = pandas.read_csv(spladderOut, compression= 'gzip', sep = '\t', dtype = str)
    as_events.rename({'contig':'chrom'}, axis = 'columns', inplace = True)

    '''
    TO DO:
    columns for mutex: exon_pre_start  exon_pre_end    exon1_start     exon1_end       exon2_start     exon2_end       exon_aft_start  exon_aft_end
    '''
    # generate columns to specify junctions between exon skips since star records intron junction not exon
    # pre_intron_start and aft_intron_end will always be the same regardless of exon1 or exon2 since they are constant
    # the only change will be the pre_intron_end and aft_intron_start for each exon
    as_events['pre_intron_start'] = as_events['exon_pre_end'].astype(int) + 1 
    as_events['aft_intron_end'] = as_events['exon_aft_start'].astype(int) - 1
    
    # variable intron
    as_events['pre_intron_end_exon1'] = as_events['exon1_start'].astype(int) - 1
    as_events['aft_intron_start_exon1'] = as_events['exon1_end'].astype(int) + 1
    as_events['pre_intron_end_exon2'] = as_events['exon2_start'].astype(int) - 1
    as_events['aft_intron_start_exon2'] = as_events['exon2_end'].astype(int) + 1
    
    # merge data frame (ignore strand since some discrpenacies with undetermined in star vs determined in spladder)
    as_events['pre_intron_start'] = as_events['pre_intron_start'].astype(str)
    as_events['aft_intron_end'] = as_events['aft_intron_end'].astype(str)
    as_events['pre_intron_end_exon1'] = as_events['pre_intron_end_exon1'].astype(str)
    as_events['aft_intron_start_exon1'] = as_events['aft_intron_start_exon1'].astype(str)
    as_events['pre_intron_end_exon2'] = as_events['pre_intron_end_exon2'].astype(str)
    as_events['aft_intron_start_exon2'] = as_events['aft_intron_start_exon2'].astype(str)
    
    # merge star intron annotation
    tmp = junctions.rename({'intron_start':'pre_intron_start', 'intron_end':'pre_intron_end_exon1', 'strandSTAR':'strandSTAR_exon1_pre', 'motif':'motif_exon1_pre', 'annotStatus':'annotStatus_exon1_pre' }, axis = 'columns')
    tmp_info_combine = as_events.merge(tmp, how = 'inner', on = ['chrom', 'pre_intron_start', 'pre_intron_end_exon1'])
    del tmp
    tmp = junctions.rename({'intron_start':'aft_intron_start_exon1', 'intron_end':'aft_intron_end', 'strandSTAR':'strandSTAR_exon1_aft', 'motif':'motif_exon1_aft', 'annotStatus':'annotStatus_exon1_aft' }, axis = 'columns')
    tmp_info_combine_1 = tmp_info_combine.merge(tmp, how = 'inner', on = ['chrom', 'aft_intron_start_exon1', 'aft_intron_end'])
    del tmp
    tmp = junctions.rename({'intron_start':'pre_intron_start', 'intron_end':'pre_intron_end_exon2', 'strandSTAR':'strandSTAR_exon2_pre', 'motif':'motif_exon2_pre', 'annotStatus':'annotStatus_exon2_pre' }, axis = 'columns')
    tmp_info_combine = tmp_info_combine_1.merge(tmp, how = 'inner', on = ['chrom', 'pre_intron_start', 'pre_intron_end_exon2'])
    del tmp
    tmp = junctions.rename({'intron_start':'aft_intron_start_exon2', 'intron_end':'aft_intron_end', 'strandSTAR':'strandSTAR_exon2_aft', 'motif':'motif_exon2_aft', 'annotStatus':'annotStatus_exon2_aft' }, axis = 'columns')
    final = tmp_info_combine.merge(tmp, how = 'inner', on = ['chrom', 'aft_intron_start_exon2', 'aft_intron_end'])
    
    # free memory
    del tmp, tmp_info_combine, tmp_info_combine_1
    gc.collect()
    
    # apply data filtering criteria
    events_of_interest = filtering_criteria(event_type = 'mutex', base_df = final, geneMatch = args.geneMatchFile, annot_filter = args.novel, diff_exp = args.deTestResults, pval_adj = args.pvalAdj)

    # only get unique event ids.-- drop all duplicated event ids and store counts of duplication in df for later annotation as warning for false positive of novel junction
    unique_events_of_interest = events_of_interest.drop(['annotStatus_exon1_pre', 'annotStatus_exon1_aft', 'annotStatus_exon2_pre', 'annotStatus_exon2_aft'], axis=1)
    unique_events_of_interest.drop_duplicates(keep = 'first', inplace = True)
    #total_events = unique_events_of_interest['event_id'].value_counts()
        
    '''
    Experimental blast implementation
    *** needs testing ***
    '''
    
    # TO DO
    
    '''
    end of experimental blast implementation
    '''  
    
    # check if strandSTAR and strand are concordant
    #       if alt1 and alt2 are not concordant within just star, then skip\
    #       if not concordant between star and strand, check if strandSTAR is ud
    #            if ud, then use strand
    #       if concordant then continue

    for idx,row in unique_events_of_interest.iterrows():
        try:
            assert row['strandSTAR_exon1_pre'] == row['strandSTAR_exon2_pre'], 'There is a strand conflict at the following event_id: {}.  Skipping this event.'.format(row['event_id'])
        except  AssertionError as err:
            print(err)
            continue
        
        
        ## TO DO -- needs reworking ##
        # args.kmer has one subtracted since the ase event is built into the upsteam and downstream contant exons    
        
        if ((row['strandSTAR_exon1_pre'] != row['strand']) & (row['strandSTAR_exon1_pre'] == 'ud')): # note it does not matter if we use strandSTAR_pre for strandStar_aft since we already checked the assertion that the STAR strands are concordant with each other
            orfs_exon_1 = calculate_orfs(event_type = args.eventType, fasta = dna_fasta, kmer_length = args.kmer, strand = row['strand'], chrom = row['chrom'], flank_left_start = int(row['exon_pre_start']), flank_left_end = int(row['exon_pre_end']), flank_right_start = int(row['exon_aft_start']), flank_right_end = int(row['exon_aft_end']), ase_start = int(row['exon1_start']), ase_end = int(row['exon1_end']))
            peptides_exon1 = translate_orfs(event_id = row['event_id'], orfs =  orfs_exon_1, peptide_bank = peptides_exon1)
            del orfs_exon_1
            gc.collect()
            orfs_exon_2 = calculate_orfs(event_type = args.eventType, fasta = dna_fasta, kmer_length = args.kmer, strand = row['strand'], chrom = row['chrom'], flank_left_start = int(row['exon_pre_start']), flank_left_end = int(row['exon_pre_end']), flank_right_start = int(row['exon_aft_start']), flank_right_end = int(row['exon_aft_end']), ase_start = int(row['exon2_start']), ase_end = int(row['exon2_end']))
            peptides_exon2 = translate_orfs(event_id = row['event_id'], orfs = orfs_exon_2, peptide_bank = peptides_exon2)
            del orfs_exon_2
            gc.collect()
            
            
        elif ((row['strandSTAR_exon1_pre'] != row['strand']) & (row['strandSTAR_exon1_pre'] != 'ud')):
            try:
                if args.testMode: # forces the translation despite mismatch or NaN and use spladder's strand since just for testing mode
                    orfs_exon_1 = calculate_orfs(event_type = args.eventType, fasta = dna_fasta, kmer_length = args.kmer, strand = row['strand'], chrom = row['chrom'], flank_left_start = int(row['exon_pre_start']), flank_left_end = int(row['exon_pre_end']), flank_right_start = int(row['exon_aft_start']), flank_right_end = int(row['exon_aft_end']), ase_start = int(row['exon1_start']), ase_end = int(row['exon1_end']))
                    peptides_exon1 = translate_orfs(event_id = row['event_id'], orfs =  orfs_exon_1, peptide_bank = peptides_exon1)
                    del orfs_exon_1
                    gc.collect()
                    orfs_exon_2 = calculate_orfs(event_type = args.eventType, fasta = dna_fasta, kmer_length = args.kmer, strand = row['strand'], chrom = row['chrom'], flank_left_start = int(row['exon_pre_start']), flank_left_end = int(row['exon_pre_end']), flank_right_start = int(row['exon_aft_start']), flank_right_end = int(row['exon_aft_end']), ase_start = int(row['exon2_start']), ase_end = int(row['exon2_end']))
                    peptides_exon2 = translate_orfs(event_id = row['event_id'], orfs = orfs_exon_2, peptide_bank = peptides_exon2)
                    del orfs_exon_2
                    gc.collect()
            except:
                    print('star strand is {} and spladder strands is {}.  Skipping event.'.format(row['strandSTAR_exon1_pre'], row['strand']))
                    continue
        
        elif (row['strandSTAR_exon1_pre'] == row['strand']):
            orfs_exon_1 = calculate_orfs(event_type = args.eventType, fasta = dna_fasta, kmer_length = args.kmer, strand = row['strandSTAR_exon1_pre'], chrom = row['chrom'], flank_left_start = int(row['exon_pre_start']), flank_left_end = int(row['exon_pre_end']), flank_right_start = int(row['exon_aft_start']), flank_right_end = int(row['exon_aft_end']), ase_start = int(row['exon1_start']), ase_end = int(row['exon1_end']))
            peptides_exon1 = translate_orfs(event_id = row['event_id'], orfs =  orfs_exon_1, peptide_bank = peptides_exon1)
            del orfs_exon_1
            gc.collect()
            orfs_exon_2 = calculate_orfs(event_type = args.eventType, fasta = dna_fasta, kmer_length = args.kmer, strand = row['strandSTAR_exon2_pre'], chrom = row['chrom'], flank_left_start = int(row['exon_pre_start']), flank_left_end = int(row['exon_pre_end']), flank_right_start = int(row['exon_aft_start']), flank_right_end = int(row['exon_aft_end']), ase_start = int(row['exon2_start']), ase_end = int(row['exon2_end']))
            peptides_exon2 = translate_orfs(event_id = row['event_id'], orfs = orfs_exon_2, peptide_bank = peptides_exon2)
            del orfs_exon_2
            gc.collect()
    
    
    # TO DO: calcuate percent overlap of kmer with flanking region
    
    
    # get psi info columns
    unique_events_of_interest.filter(regex = '.psi', axis = 1)
    # get read coverage info columns
    unique_events_of_interest.filter(regex = '_cov', axis = 1)
    # get mean event count and mean gene expression columns
    unique_events_of_interest.filter(regex = '^mean_', axis = 1)
    # get log2FC of differential testing in event counts and gene expression levels columns
    unique_events_of_interest.filter(regex = '^log2FC_', axis = 1)
        

def multiExon_skip():
    '''
    input:
    description:
    return:
    '''
    # should be the same as exon_skip
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
    parser.add_argument('--testMode', action = 'store_true', help = 'Ignores some of the filtering criteria so all input cases are used for testing code.')
    parser.add_argument('--knownJunctions', type = str, required = True, help = 'path to sjdbList.fromGTF.out.tab generated within STARgenome directory' )
    parser.add_argument('--geneMatchFile', type = str, help = 'The output of greping the 9th column of a gtf file (see extras folder to find grep command) or a tab-delimited ensembl to gene name match')
    parser.add_argument('--blast', type = str, help = 'Full path to executable, including executable to standalone blastp software from NCBI')
    parser.add_argument('--blastDb', type = str, help = 'Full path to blastp database to use for identifying proper protein reading frame; typically this is a genome-wide uniport annotated protein database formatted for standalone BLAST for your organism; see extras folder for information on how to generate this')
    parser.add_argument('--frameMatch', action = 'store_true', help = 'only translates and returns peptide that have a curated protein match from blast results in correct reading frame.  Otherwise, returns all frames regardless of match.')
    
    args = parser.parse_args()
    
    codon_library = generate_codon_reference()
    dna_fasta = SeqIO.index(args.fasta, 'fasta')
    junctions = combine_data(junctionFiles = args.SJfiles, annotation = args.knownJunctions)
    
    if args.eventType == 'intron_retention':
        intron_retention(junctions = junctions, spladderOut = args.ase, outdir = args.outdir)
        
    if args.eventType == 'exon_skip':
        exon_skip(junctions = junctions, spladderOut = args.ase, outdir = args.outdir)

    if args.eventType == 'three_prime_alt':
        three_prime_alt(junctions = junctions, spladderOut = args.ase, outdir = args.outdir)
        
    if args.eventType == 'five_prime_alt':
        five_prime_alt(junctions = junctions, spladderOut = args.ase, outdir = args.outdir)
    
    if args.eventType == 'mutex':
        mutex(junctions = junctions, spladderOut = args.ase, outdir = args.outdir)
