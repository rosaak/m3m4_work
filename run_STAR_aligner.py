#!/usr/bin/env python3
# -*- coding: utf-8 -*-

__author__ ="roshan padmanbhan"
__version__=0.1

"""
Written for m3m4 mouse model rnaseq

"""

import os
from  pathlib import Path
import shlex
import subprocess
from collections import OrderedDict
import logging
import time


def _run_star_aligner( threads,genome_dir, out_prefix, R1, R2, zip=True ):
    zip_option = ' --readFilesCommand gunzip -c '
    if zip :
        cmd = shlex.split('STAR --runThreadN ' + str(threads) + ' --genomeDir ' + genome_dir + ' --outFileNamePrefix ' + out_prefix + zip_option + ' --readFilesIn ' + R1 + ' ' + R2  )
    else :
        cmd = shlex.split('STAR --runThreadN ' + str(threads) + ' --genomeDir ' + genome_dir + ' --outFileNamePrefix ' + out_prefix + ' --readFilesIn ' + R1 + ' ' + R2 )
    return(cmd)

def _make_logger( namex ):
    logger = logging.getLogger(__name__)
    logger.setLevel(logging.INFO)
    # create a file handler
    handler = logging.FileHandler('star_'+namex+'_'+time.ctime().replace(' ','_').replace(':','-')+'.log')
    handler.setLevel(logging.INFO)
    formatter = logging.Formatter('%(asctime)s - %(message)s')
    handler.setFormatter(formatter)
    logger.addHandler(handler)
    return( logger )

def _make_dir( path_x ):
    td = Path( path_x )
    if not td.is_dir():
        td.mkdir()
    return( path_x )

def _run_commands( cmd ):
    s = subprocess.call(cmd, shell=False)
    return(s)

def run_star(tissue, dictx, genome_index, res_loc):
    log = _make_logger( tissue )
    log.info('Start STAR Alignment\n')
    res_loc = _make_dir( os.path.join( res_loc, tissue ) )
    log.info( 'Result Location dir {}\n'.format( res_loc ))
    for i in dictx.keys():
        R1 = dictx.get( i )
        R2 = str(dictx.get( i )).replace('_R1','_R2')
        bn = dictx.get( i ).name
        bn = bn.split('.')[4].split('_')[0]
        prefix = os.path.join( res_loc, bn )
        log.info('Processing files {}'.format(bn) )
        log.info('Result prefix {}'.format(prefix) )
        if R1.name.endswith('.gz'):
            log.info( ' '.join( _run_star_aligner( threads=15, genome_dir=genome_index, out_prefix=prefix,  R1=str(R1), R2=R2, zip=True ))+"\n" )
            _run_commands(_run_star_aligner( threads=15, genome_dir=genome_index, out_prefix=prefix,  R1=str(R1), R2=R2, zip=True  ))
        else :
            log.info( ' '.join( _run_star_aligner( threads=15, genome_dir=genome_index, out_prefix=prefix,  R1=str(R1), R2=R2, zip=False ))+"\n")
            _run_commands(_run_star_aligner( threads=15, genome_dir=genome_index, out_prefix=prefix,  R1=str(R1), R2=R2, zip=False ))
    log.info( 'STAR aligner finished\n' )

def main():
    genome_star_index_fp='/home/padmanr/niazif-share/Library_Files/StarIndex/STARindex/mm10'
    fastq_loc='/mnt/isilon/data/w_gmi/gmi-to-be-archived/engclab_ngs/m3m4_p14brain_p40cortex_p40thyroid_RNA-Seq_2013/Raw_fastq'
    res_loc_fp = '/home/padmanr/niazif-share/Stetson/Amanda_M3M4_Reanalysis/data/star_aligned'

    # make some dicts
    files_dict = OrderedDict()
    fqx = Path(fastq_loc)
    files_dict = { i.name : i for i in sorted( fqx.glob('*_R1.*')) }

    thyroid_dict = { i : files_dict.get(i) for i in files_dict.keys() if 'thyroid' in i  }
    brain_dict = { i : files_dict.get(i) for i in files_dict.keys() if 'brain' in i  }
    cortex_dict = { i : files_dict.get(i) for i in files_dict.keys() if 'cortex' in i  }

    run_star( tissue = 'cortex', dictx=cortex_dict, genome_index=genome_star_index_fp, res_loc=res_loc_fp )
    run_star( tissue = 'brain', dictx=brain_dict, genome_index=genome_star_index_fp, res_loc=res_loc_fp )
    run_star( tissue = 'thyroid', dictx=thyroid_dict, genome_index=genome_star_index_fp, res_loc=res_loc_fp )

if __name__ == '__main__':
    main()
