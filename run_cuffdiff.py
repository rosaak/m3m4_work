#!/usr/bin/env python3
# -*- coding: utf-8 -*-

__author__ ="roshan padmanbhan"
__version__=0.1

"""
Written for m3m4 mouse model rnaseq
Run cuffdiff on the combinations
"""

import os
from  pathlib import Path
import shlex
import subprocess
import logging
import time
import sys

"""
Aim : run cuffdiff
    cuffdiff \
    -o diff_out \
    -b genome.fa \
    -p 8 \
    –L C1,C2 \
    -u merged_asm/merged.gtf \
    ./C1_R1_thout/accepted_hits.bam,./C1_R2_thout/accepted_hits.bam,./C1_R3_thout/accepted_hits.bam \
    ./C2_R1_thout/accepted_hits.bam,./C2_R3_thout/accepted_hits.bam,./C2_R2_thout/accepted_hits.bam

genome
    mm10
labels
    1. het_mut
    2. het_wt
    3. wt_mut
merged files for
    1. het , mut 
    2. het , wt
    3. wt , mut
"""

# all functions
def run_commands( cmd ):
    s = subprocess.call(cmd, shell=False)
    return(s)

def _cuffdiff_v1( labels, threads, outdir, genome, mergedgtf, accepted_hits ):
    """
    labels = comma sep string of labels ex: Het,Mut
    accepted_hits = str of bam files from condition A and condition B
    """
    cmd = shlex.split("cuffdiff -L " + labels + " -p " + str(threads) + " -o "+ outdir + "  -b " + genome + "  -u " + mergedgtf + " " + accepted_hits )
    return(cmd)

def get_accepted_hits( day, L1, L2, l1, l2):
    """str,str,str,list,list -> str
    L1, L2 : label for conditions
    l1 and l2 : list location of conditions
    Returns a string of locations to the accepted_hits.bam files  
    """
    assemblies = l1 + l2
    L1x=[]
    L2x=[]
    for i in assemblies:
        if day in i:
            if L1 in str(i):
                L1x.append( str(i) )
            if L2 in str(i):
                L2x.append( str(i) )
    L1x = ','.join( L1x )
    L2x = ','.join( L2x )
    return( L1x+' '+ L2x)

def make_logger( ):
    logger = logging.getLogger(__name__)
    logger.setLevel(logging.INFO)
    # create a file handler
    handler = logging.FileHandler('cuffdiff'+time.ctime().replace(' ','_').replace(':','-')+'.log')
    handler.setLevel(logging.INFO)
    formatter = logging.Formatter('%(asctime)s - %(message)s')
    handler.setFormatter(formatter)
    logger.addHandler(handler)
    return( logger )

log = make_logger( )
log.info('Start cuffdiff')

# locations
cufflinks_loc="/home/padmanr/niazif-share/Stetson/Amanda_M3M4_Reanalysis/tophat/cufflinks"
bam_loc="/home/padmanr/niazif-share/Marilyn/Stetson_bam"
cuffmerged_loc="/home/padmanr/niazif-share/Stetson/Amanda_M3M4_Reanalysis/tophat/cuffmerge_03"
cuffdiff_res_loc="/home/padmanr/niazif-share/Stetson/Amanda_M3M4_Reanalysis/tophat/cuffdiff_01"
genome_mm10="/home/padmanr/niazif-share/Farshad/LIBRARY_FILES/mm10.fa"

# get all the acceptedhits bam files into a dict
clocx = Path(bam_loc)
if clocx.is_dir():
    bams  = { str(i.name).split('_')[0] : i for i in  sorted(clocx.glob('P*.bam')) }
    P14_Het_loc = [ str(i) for i in  sorted(clocx.glob( 'P*Het*.bam' ) ) if 'P14' in str(i) ]
    P14_Wt_loc  = [ str(i) for i in  sorted(clocx.glob( 'P*Wt*.bam' ) ) if 'P14' in str(i) ]
    P14_Mut_loc = [ str(i) for i in  sorted(clocx.glob( 'P*Mut*.bam' ) ) if 'P14' in str(i) ]
    P40_Het_loc = [ str(i) for i in  sorted(clocx.glob( 'P*Het*.bam' ) ) if 'P40' in str(i) ]
    P40_Wt_loc  = [ str(i) for i in  sorted(clocx.glob( 'P*Wt*.bam' ) ) if 'P40' in str(i) ]
    P40_Mut_loc = [ str(i) for i in  sorted(clocx.glob( 'P*Mut*.bam' ) ) if 'P40' in str(i) ]
else:
    sys.exit("abort dir invalid\n")

Loc_dict = {'P14_Het':P14_Het_loc, 'P14_Wt':P14_Wt_loc, 'P14_Mut':P14_Mut_loc, 'P40_Het':P40_Het_loc, 'P40_Wt':P40_Wt_loc, 'P40_Mut':P40_Mut_loc }

# for each merged gtf file run the cuffdiff
merged_gtfs = Path( cuffmerged_loc).glob('P*/merged.gtf')
for i in sorted( merged_gtfs ):
    log.info('processing {}'.format( str(i)))
    day, L1, L2 = str(i.parents[0]).split( '/' )[-1].replace('_merged','').split('_')
    c1 = [ Loc_dict.get(s) for s in Loc_dict.keys() if day+'_'+L1 in s ][0]
    c2 = [ Loc_dict.get(s) for s in Loc_dict.keys() if day+'_'+L2 in s ][0]
    labelx = L1 + ',' + L2
    out_dir = os.path.join( cuffdiff_res_loc , day+'_'+L1+'_'+L2+'_cuffdiff' )
    ahits_str = get_accepted_hits( day=day, L1=L1, L2=L2, l1=c1, l2=c2)
    log.info(' '.join( _cuffdiff_v1( labels=labelx, threads=15, outdir=out_dir, genome=genome_mm10, mergedgtf=str(i), accepted_hits=ahits_str)))
    log.info(' ')
    run_commands( _cuffdiff_v1( labels=labelx, threads=15, outdir=out_dir, genome=genome_mm10, mergedgtf=str(i), accepted_hits=ahits_str) )

log.info("finsihed running cuffdiff\n")
