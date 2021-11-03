import Bio
import pandas
import argparse

def combine_data(junctionFiles):
    
    # setting args.SJfiles is tmp just so we don't have to use parser right away for development and testing"
    junctionFiles = ['/mnt/IM_drive/Jill_Slansky/neoantigen_alternative_splicing_project_07122021/new_analysis_tonya_07212021/RNA_differential_alternative_splicing_spladder/output/1_EO771_EV_untreated_S7_spladder_inputSJ.out.tab', 
                     '/mnt/IM_drive/Jill_Slansky/neoantigen_alternative_splicing_project_07122021/new_analysis_tonya_07212021/RNA_differential_alternative_splicing_spladder/output/2_EO771_EV_untreated_S8_spladder_inputSJ.out.tab', 
                     '/mnt/IM_drive/Jill_Slansky/neoantigen_alternative_splicing_project_07122021/new_analysis_tonya_07212021/RNA_differential_alternative_splicing_spladder/output/3_EO771_EV_untreated_S9_spladder_inputSJ.out.tab', 
                     '/mnt/IM_drive/Jill_Slansky/neoantigen_alternative_splicing_project_07122021/new_analysis_tonya_07212021/RNA_differential_alternative_splicing_spladder/output/5_EO771_EV_plus_doxy_S10_spladder_inputSJ.out.tab', 
                     '/mnt/IM_drive/Jill_Slansky/neoantigen_alternative_splicing_project_07122021/new_analysis_tonya_07212021/RNA_differential_alternative_splicing_spladder/output/6_EO771_EV_plus_doxy_S11_spladder_inputSJ.out.tab', 
                     '/mnt/IM_drive/Jill_Slansky/neoantigen_alternative_splicing_project_07122021/new_analysis_tonya_07212021/RNA_differential_alternative_splicing_spladder/output/7_EO771_miR200c_plus_doxy_S12_spladder_inputSJ.out.tab', 
                     '/mnt/IM_drive/Jill_Slansky/neoantigen_alternative_splicing_project_07122021/new_analysis_tonya_07212021/RNA_differential_alternative_splicing_spladder/output/8_EO771_miR200c_plus_doxy_S13_spladder_inputSJ.out.tab', 
                    '/mnt/IM_drive/Jill_Slansky/neoantigen_alternative_splicing_project_07122021/new_analysis_tonya_07212021/RNA_differential_alternative_splicing_spladder/output/9_EO771_miR200c_plus_doxy_S14_spladder_inputSJ.out.tab'
                    ]
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



def intron_retention(junctions, spladderOut):
    spladderOut = '/mnt/IM_drive/Jill_Slansky/neoantigen_alternative_splicing_project_07122021/new_analysis_tonya_07212021/RNA_differential_alternative_splicing_spladder/spladder_out/merge_graphs_intron_retention_C2.confirmed.txt.gz'
    
    # read in data from spladder and then rename columns for df merge
    as_events = pandas.read_csv(spladderOut, compression= 'gzip', sep = '\t', dtype = str)
    as_events.rename({'contig':'chrom'}, axis = 'columns', inplace = True)
    
    # merge data frame (ignore strand since some discrpenacies with undetermined in star vs determined in spladder)
    info_combine = as_events.merge(junctions, how = 'left', on = ['chrom', 'intron_start', 'intron_end'])
    info_combine.to_csv("testing_introns.tsv", sep = "\t", index = False)
    
    #TO DO: check strand
   

def exon_skip():
    pass

def three_prime_alt():
    pass

def five_prime_alt():
    pass

def mutex():
    pass

def multiExon_skip():
    pass

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description= '')
    
    parser.add_argument('--SJfiles', required = True, action = 'extend', nargs= '+', type = str, help = 'list of all SJ.out.tab files from STAR two-pass mode')
    
    args = parser.parse_args()
    
    junctions = combine_data(junctionFiles = args.SJfiles)