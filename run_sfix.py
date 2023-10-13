#!python3

import numpy as np
import argparse, time, os, pickle, datetime
import StringFix.core as sf
from StringFix.core  import StringFix_analysis, StringFix_synthesis, Generate_Reference, StringFix_GFFRead

## file_name_wo_ext = StringFix_analysis(sam_file, gtf = None, suffix = 'stringfix', jump = 0, sa = 2, cbyc = False, \
##                    n_cores = 1, min_frags_per_gene = 2, len_th = 200, out_tr = False )
## file_name_wo_ext, file_name_genome = StringFix_synthesis(gtf_file, genome_file, rev = True, n_cores = 1)
## file_name_ref_tr, file_name_ref_pr = Generate_Reference(gtf_file, genome_file)
## file_name_wo_ext = StringFix_GFFRead(gtf_file, genome_file = None, n_cores = 1)
## df, dft, df_sel = run_blast( inpt_fa, ref_fa, path_to_blast = None, trareco = True, ref_info = False, \
##            dbtype = 'nucl', ref = 0, sdiv = 10, mx_ncand = 6, verbose = False)

def get_args():
    parser = argparse.ArgumentParser(description='', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser._action_groups.pop()
    required = parser.add_argument_group('REQUIRED PARAMETERS')
    optional = parser.add_argument_group('OPTIONAL PARAMETERS')

    required.add_argument('-input', type=str, metavar='',
                          help='Input SAM or BAM file. pybam.py is required for BAM input.')
    optional.add_argument('-gnm', type=str, metavar='',
                          help='Genome fasta file used to generate the input SAM/BAM. It is used to correct ' + \
                            'sequences not covered by reads, to identify SNPs and to generate customized genome.')
    optional.add_argument('-gtf', type=str, metavar='',
                          help='GTF/GFF file to be used for guide annotation ')
    optional.add_argument('-dir', type=str,  metavar='',
                        help='Output directory', default='SFix_out')
    optional.add_argument('-p', type=str,  metavar='',
                        help='Number of cores to use', default="4")
    optional.add_argument('-g', action = 'store_true',
                        help='Use this option to generate customized genome.')
    # optional.add_argument('-c', action = 'store_true',
    #                     help='Set this option to process chromosomes 1-by-1 to save memory usage.')
    optional.add_argument('-sfx', type=str,  metavar='',
                        help='Suffix to be used for output file names.', default = None)
    optional.add_argument('-mcv', type=str,  metavar='',
                        help='Minimum read coverage for inclusion in the output GTF/GFF and transcriptome/proteome [0~1]. For a length L transcript, \
                              it will be included in the output if at least mcv*L bases are covered by reads', default="0.5")
    optional.add_argument('-mcd', type=str,  metavar='',
                        help='Minimum base coverage depth for inclusion in the output', default="1")
    optional.add_argument('-mtl', type=str,  metavar='',
                        help='Minimum transcript length in bases.', default="200")
    optional.add_argument('-mpl', type=str,  metavar='',
                        help='Minimum protein (amino acid sequence) length.', default="50")
    optional.add_argument('-mdi', type=str,  metavar='',
                        help='Minimum coverage depth for an insersion to be considered valid [>=1.0].', default="2")
    optional.add_argument('-mdd', type=str,  metavar='',
                        help='Minimum coverage depth for a deletion to be considered valid [>=1.0].', default="2")
    optional.add_argument('-mdp', type=str,  metavar='',
                        help='Minimum coverage depth for SNP correction [>=1.0].', default="3")
    optional.add_argument('-mdf', type=str,  metavar='',
                        help='Minimum coverage depth fraction for insersion and deletion to be considered valid (0~1).', default="0.25")
    optional.add_argument('-np', type=str,  metavar='',
                        help='Maximum number of isoforms per gene(read chunk). Applies when guide GTF is not provided. (>8).', default="16")
    optional.add_argument('-xsa', action = 'store_true',
                        help='Use this option to exclude secondary alignment.')
    optional.add_argument('-s', action = 'store_false',
                        help='Use this option to save gene descriptor and loci.', default=False)
    optional.add_argument('-j', type=str,  metavar='',
                        help='Jump to step j [0 or 1 or 2]', default="0")
    args = parser.parse_args()
    return args, parser


def main():

    print('+------------------------------+')
    print('|      Stringfix v.0.7.1       |')
    print('+------------------------------+')

    args, parser = get_args()

    if (args.input is None): # | (args.gtf is None) | (args.gnm is None):
        parser.print_help()
        return
    else:

        genome_fa = args.gnm
        sam_file = args.input
        gtf_file = args.gtf
        n_cores = int(args.p)

        if (genome_fa is not None) & (gtf_file is not None):
            print('StringFix is running in annotation guided mode with sequence correction.')
        elif (genome_fa is not None) & (gtf_file is None):
            print('StringFix is running in genome guided mode with sequence correction.')
        elif (genome_fa is None) & (gtf_file is None):
            print('StringFix is running in genome guided mode without sequence correction.')
        else:
            print('StringFix is running in annotation guided mode without sequence correction.')
        

        # rev = False
        rev = args.g
        cbyc = True # args.c
        sa = 2
        if args.xsa: sa = 0

        save_gds_and_loci = False
        if args.s: save_gds_and_loci = True

        sf.MIN_CVG_FOR_SNP_CORRECTION = int(args.mdp)
        sf.I_VALID_CVG = max(1, float(args.mdi))
        sf.D_VALID_CVG = max(1, float(args.mdd))
        sf.I_VALID_CFRAC = float(args.mdf)
        sf.D_VALID_CFRAC = float(args.mdf)
        sf.MIN_COV_TO_SEL = float(args.mcv)
        sf.MIN_ABN_TO_SEL = float(args.mcd)
        sf.MIN_TR_LENGTH = int(args.mtl)
        sf.MIN_PR_LENGTH = int(args.mpl)
        sf.MAX_NUM_PATHS = int(args.np)
        sf.GEN_AAS_UNSPEC = True
        jump_to = int(args.j)

        if gtf_file is None:
            sf.MAX_NUM_PATHS = 20
            sf.MIN_N_RDS_PER_GENE = 30
            sf.TR_FILTER_SEL = 2

        if args.dir is not None: 
            if not os.path.exists(args.dir): os.mkdir(args.dir)

        file_out, gds = StringFix_analysis(sam_file, gtf = gtf_file, genome_file = genome_fa, jump = jump_to, \
                n_cores=n_cores, out_tr = True, cbyc = cbyc, suffix = args.sfx, out_dir = args.dir, sa = sa, \
	  sav_rgn =  save_gds_and_loci) 

        if (save_gds_and_loci):  
            print('Saving gene descriptor .. ', end = ' ', flush = True)
            with open(file_out + '.gds', 'wb') as f:
                pickle.dump( (gds), f )
            print('done. ')
        elif (gds is None) | (jump_to >= 2):
            print('StringFix started ', datetime.datetime.now())
            print('Loading gene descriptor .. ', end = ' ', flush = True)
            with open(file_out + '.gds', 'rb') as f:
                gds = pickle.load( f )
            print('done. ')
        
        if (genome_fa is None):
            print('StringFix finished ', datetime.datetime.now())
        else:
            file_out, file_out_genome = StringFix_synthesis( file_out + '.gff', genome_fa, \
                    n_cores=n_cores, rev = rev, out_dir = args.dir, gds = gds )

if __name__=="__main__":
    main()

