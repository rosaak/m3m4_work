#!/usr/bin/env python3
# -*- coding: utf-8 -*-

__author__ ="roshan padmanabhan <rosaak@gmail.com>"
__version__=0.1

"""

Written for m3m4 mouse model rnaseq

Aim : run cuffmerge 
merge the cleanup cufflinks output
there will be three mergings
1. het , mut 
2. het , wt
3. wt , mut
"""

import os
from  pathlib import Path
import shlex
import subprocess
from collections import OrderedDict
import sys

def _make_dir( path_x ):
    td = Path( path_x )
    if not td.is_dir():
        td.mkdir()

def run_commands( cmd ):
    s = subprocess.call(cmd, shell=False)
    return(s)

def _cuffmerge( gtf, threads, outdir, assembly_txt ):
    cmd = shlex.split("cuffmerge -g " + gtf + " -p " + str(threads) + " -o "+ outdir + " " + assembly_txt)
    return(cmd)

def _cuffmerge_v2( threads, outdir, assembly_txt ):
    cmd = shlex.split("cuffmerge --keep-tmp -p " + str(threads) + " -o "+ outdir + " " + assembly_txt)
    return(cmd)

def _cuffmerge_v3( gtf, threads, outdir, genome, assembly_txt ):
    cmd = shlex.split("cuffmerge --keep-tmp -g " + gtf + " -p " + str(threads) + " -o "+ outdir + "  -s " + genome + " " + assembly_txt )
    return(cmd)

def write_assembly_file( day, aname , c1 , c2, transcripts ):
    assemblies = c1 +  c2
    assembly_txt_fp = os.path.join( assembly_loc, aname +'_assembly.txt')
    to_write = list()
    for i in assemblies :
        bn = i.split('/')[-1]
        if day in bn:
            to_write.append( str(transcripts.get( bn ))+'\n' )
        with open( assembly_txt_fp, 'w') as ofh :
                ofh.writelines( to_write )
    return( assembly_txt_fp )

# locations from hpc 
loc="/home/padmanr/niazif-share/Stetson/Amanda_M3M4_Reanalysis/tophat/cufflinks"
assembly_loc="/home/padmanr/niazif-share/Stetson/Amanda_M3M4_Reanalysis/tophat/cuffmerge_04"
mm10_gtf="/home/padmanr/niazif-share/Roshan/bank/Mus_musculus/UCSC/mm10/Annotation/Genes/genes.gtf"
genome_mm10="/home/padmanr/niazif-share/Farshad/LIBRARY_FILES/mm10.fa"

# get all the transcript file loations
locx=Path(loc)
if locx.is_dir():
    transcripts = OrderedDict()
    transcripts = { str(i.parents[0]).split('/')[-1] : i for i in  sorted(locx.glob('P*/transcripts.gtf')) }
else:
    sys.exit("abort dir invalid\n")

# Make assembly.txt files and save them
_make_dir(assembly_loc)

# get all het , wt abd mut 
Het_loc = [ str(i) for i in  sorted(locx.glob( 'P*Het*') )]
Wt_loc  = [ str(i) for i in  sorted(locx.glob( 'P*Wt*') )]
Mut_loc = [ str(i) for i in  sorted(locx.glob( 'P*Mut*') )]

# making three cuffmerges assembly.txt files, all het_mut, het_wt and wt_mut for day 14 and day 40
assembly_list=list()
day = "P14"
assembly_list.append( write_assembly_file(day,'P14_Het_Mut', Het_loc, Mut_loc, transcripts) )
assembly_list.append( write_assembly_file(day,'P14_Het_Wt', Het_loc, Wt_loc, transcripts) )
assembly_list.append( write_assembly_file(day,'P14_Wt_Mut', Wt_loc, Mut_loc, transcripts) )

day = "P40"
assembly_list.append( write_assembly_file(day,'P40_Het_Mut', Het_loc, Mut_loc, transcripts) )
assembly_list.append( write_assembly_file(day,'P40_Het_Wt', Het_loc, Wt_loc, transcripts) )
assembly_list.append( write_assembly_file(day,'P40_Wt_Mut', Wt_loc, Mut_loc, transcripts) )

# run cuffmerge 
for i in assembly_list:
    print( "processing : ", i.split("/")[-1].replace('_assembly.txt',''), end="\n" )
    res_loc = os.path.join( assembly_loc , i.split("/")[-1].replace('_assembly.txt','_merged') )
    print( ' '.join( _cuffmerge_v3( gtf=mm10_gtf, genome=genome_mm10, threads=15, outdir=res_loc, assembly_txt=i ) ))
    print()
    run_commands( _cuffmerge_v3( gtf=mm10_gtf, genome=genome_mm10, threads=15, outdir=res_loc, assembly_txt=i ) )

print("done.....")
