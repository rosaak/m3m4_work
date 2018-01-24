from pathlib import Path
import os
import pandas as pd
import shutil

"""
Date : 28, Dec, 2017
This script takes the cufflinks res dir and go tru each file and removes non cannonical chr files
I moved the skipped.gtf to new locations
"""

def clean_transcripts( file_transcript ):
    """
    return a new df with canonical df
    """
    header=['seqname','source','feature','start','end','score','strand','frame','attributes']
    transcripts =  pd.read_csv( file_transcript , sep="\t", names=header)
    can_chr = sorted( [ i for i in transcripts.seqname.value_counts().index.tolist() if 'chr' in str(i) and '_' not in str(i) ]  )
    df_filtered = [ ]
    for i in can_chr:
        tmp_df = transcripts[ transcripts.seqname == i]
        df_filtered.append( tmp_df )
    df_new = pd.concat( df_filtered )
    return( df_new )

def clean_fpkm_tracking( file_transcript ):
    """
    return a new df with canonical df
    """
    f_tracking =  pd.read_csv( file_transcript , sep="\t")
    can_chr = sorted( [ i for i in pd.DataFrame(f_tracking.locus.str.split(':',1).tolist(), columns =["chr", "loc"] ).chr.value_counts().index.tolist() if 'chr' in str(i)  and '_' not in str(i) ])
    df_filtered = [ ]
    for i in can_chr:
        tmp_df = f_tracking[ pd.DataFrame(f_tracking.locus.str.split(':',1).tolist(), columns =["chr", "loc"] ).chr == i ]
        df_filtered.append( tmp_df )
    df_new = pd.concat( df_filtered )
    return( df_new )

def save_df( df, res_fp, sep, header ):
    """
    write / overwrites the pandas df to res_fp
    """
    res_fpx = Path( res_fp )
    if res_fpx.is_file():
        df.to_csv( str( res_fpx ), index=False, index_label=False,  header=header, sep=sep, quoting=0 , quotechar=" ")
    else :
        res_fpx.touch()
        df.to_csv( str( res_fpx ) , index=False, index_label=False,  header=header, sep=sep, quoting=0 , quotechar=" ")
    if res_fpx.stat().st_size >= 0 :
        return( "df saved to {}".format( str( res_fpx ) ))

def main():
    # run through all cufflinks_output2 dir
    # get the file and clean up and save then in a new dir
    res_dir_fp="/home/padmanr/Nia/Stetson/Amanda_M3M4_Reanalysis/tophat/cufflinks_output_cleaned"
    input_dir_fp="/home/padmanr/Nia/Stetson/Amanda_M3M4_Reanalysis/tophat/cufflinks_output2/"
    input_dir_fpx =  Path( input_dir_fp )
    for i in sorted(input_dir_fpx.glob('*')):
        locx = Path( i )
        print('Processing : {}'.format( str(locx)) )
        indiv_res_dir = os.path.join( res_dir_fp, str( locx.name) )
        indiv_res_dirx = Path(indiv_res_dir)
        indiv_res_dirx.mkdir()
        # res locations
        transcript_res_fp = os.path.join( indiv_res_dir , "transcripts.gtf")
        isoform_res_fp = os.path.join( indiv_res_dir , "isoforms.fpkm_tracking")
        genes_res_fp = os.path.join( indiv_res_dir , "genes.fpkm_tracking")
        skipped_res_fp = os.path.join( indiv_res_dir , "skipped.gtf")
        # input locations
        transcript_loc =  [ str(i) for i in locx.glob('transcripts.gtf') ][0]
        isoformfpkm_loc =  [ str(i) for i in locx.glob('isoforms.fpkm_tracking') ][0] 
        genesfpkm_loc = [ str(i) for i in locx.glob('genes.fpkm_tracking') ][0]
        skipped_loc = [ str(i) for i in locx.glob('skipped.gtf') ][0]
        # clean up 
        save_df(df=clean_transcripts( transcript_loc ), res_fp=transcript_res_fp, sep="\t", header=False ) 
        save_df(df=clean_fpkm_tracking( isoformfpkm_loc ), res_fp=isoform_res_fp, sep="\t", header=True ) 
        save_df(df=clean_fpkm_tracking( genesfpkm_loc ), res_fp=genes_res_fp, sep="\t", header=True )
        shutil.copy( skipped_loc, skipped_res_fp)

if __name__ == '__main__':
    main()
