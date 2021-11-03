import Bio
import pandas
import argparse

def combine_data(junctionFiles):
    
    # setting args.SJfiles is tmp just so we don't have to use parser right away for development and testing"
    args.SJfiles = ['/mnt/IM_drive/Jill_Slansky/neoantigen_alternative_splicing_project_07122021/new_analysis_tonya_07212021/RNA_differential_alternative_splicing_spladder/output/1_EO771_EV_untreated_S7_spladder_inputSJ.out.tab', 
                    '/mnt/IM_drive/Jill_Slansky/neoantigen_alternative_splicing_project_07122021/new_analysis_tonya_07212021/RNA_differential_alternative_splicing_spladder/output/5_EO771_EV_plus_doxy_S10_spladder_inputSJ.out.tab', 
                    '/mnt/IM_drive/Jill_Slansky/neoantigen_alternative_splicing_project_07122021/new_analysis_tonya_07212021/RNA_differential_alternative_splicing_spladder/output/9_EO771_miR200c_plus_doxy_S14_spladder_inputSJ.out.tab'
                    ]
    # read in all SJ.tab files detected by 2-pass STAR and get the columns: chrom, intron start, intron end, strand, novelty
    junctionDFs = [pandas.read_table(inputFile, sep = '\t', header = None, names = ['chrom', 'intronStart', 'intronEnd', 'strand', 'motif', 'annotStatus'], 
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
    junctions.replace({'strand': {0:'ud' , 1: '+', 2:'-'}, 
                       'motif':{0: 'nc', 1:'GT/AG', 2:'CT/AC', 3:'GC/AG', 4:'CT/GC', 5:'AT/AC', 6:'GT/AT'},
                       'annotStatus':{0:'novel', 1:'annotated'}}, 
                       inplace = True)
    
    return junctions


def intron_retention():
    pass

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