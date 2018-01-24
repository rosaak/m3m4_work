# Aim : cuffdiff to cumbersplicR plots
# using expression cut off : 10
# using only these filters c('geneOK', 'expressedGenes', 'sigGenes', 'isoOK', 'expressedIso', 'isoClass', 'sigIso' )
# version 01


suppressPackageStartupMessages( library( spliceR )  )
require("BSgenome.Mmusculus.UCSC.mm10",character.only=TRUE)

msg_print = function( log, msg ){
    print( paste0(log,": ",as.character( msg )) )
}

make_cummeRbund_plot_fp = function( splr_loc_x, fn, ext ){
    return( paste0( splr_loc_x,'/',fn, ext ) )
}

make_spliceR_plot_fp = function( splr_loc_x, fn, i, x,  ext ){
    return( paste0( splr_loc_x,'/',fn, "_", i, "_", x ,ext ) )
}

make_fp = function( splr_loc_x, fn ){
    return( paste0( splr_loc_x,'/',fn) )
}

run_spliceR_plots = function(slist, evaluate, asType, out_pdf_fp ){
    pdf(out_pdf_fp)
    xx = spliceRPlot(slist, evaluate=evaluate, asType=asType)
    dev.off()
    return(xx)
}

# get all the cuffdiff location
cuffdiff_loc='/home/padmanr/niazif-share/Stetson/Amanda_M3M4_Reanalysis/tophat/cuffdiff_01'
cuffdiff_dirs = list.dirs( cuffdiff_loc,recursive=F)
cuffdiff_dir_names = list.dirs( cuffdiff_loc ,recursive=F, full.names=F)
names( cuffdiff_dirs ) = cuffdiff_dir_names

# merged gtf locations
cuffmerge_loc='/home/padmanr/niazif-share/Stetson/Amanda_M3M4_Reanalysis/tophat/cuffmerge_03'
cuffmerge_dirs = list.dirs( cuffmerge_loc,recursive=F)
cuffmerge_dir_names = list.dirs( cuffmerge_loc ,recursive=F, full.names=F)
names( cuffmerge_dirs ) = cuffmerge_dir_names

# Copy the files to a new dir
res_loc = '/home/padmanr/niazif-share/Stetson/Amanda_M3M4_Reanalysis/tophat/splicr_exp_100_jan24/' # edit here
dir.create( res_loc )
new_dirs =  sapply(cuffdiff_dir_names, function(x){ paste0( res_loc, x) })
# make new dirs
sapply( new_dirs, function(x){ dir.create(x) } )

# dfx is a data frame containing cuffdiff and new dir
dfx = as.data.frame( cbind( cuffdiff_dirs, new_dirs, cuffmerge_dirs) )
# for each cuffdiff_dir_names search in dfx then copy the files to new locations
lapply( cuffdiff_dir_names, function(x){ file.copy( list.files( as.vector( dfx[x,][,1] ), full.names=T), as.vector( dfx[x,][,2] ) )})

# mm10 annotation for making gtf files
mm10_ucscCDS = getCDS( selectedGenome="mm10", repoName="UCSC")

# run the combination of these 
expressionCutoff = 1000 # change here
evalv = c('nr_transcript','nr_genes','nr_transcript_pr_gene', 'nr_AS', 'mean_AS_gene', 'mean_AS_transcript', 'mean_transcript_exp','mean_gene_exp')
visType = c('ESI','MEE','MESI','ISI','A5','A3','ATSS','ATTS','All')
#filters = c('geneOK', 'expressedGenes', 'sigGenes', 'isoOK', 'expressedIso', 'isoClass', 'sigIso' )
filters = c('expressedGenes','geneOK', 'isoOK', 'expressedIso', 'isoClass')

# Run tru each rows of dfx
for( i in 1:length( rownames( dfx ))){ 
    nX = as.character( rownames(dfx)[i])
    nXb = gsub('_cuffdiff','', nX)
    Cond_1 = sapply(strsplit(nX ,'_'),'[[',2)
    Cond_2 = sapply(strsplit(nX ,'_'),'[[',3)
    splr_loc = as.vector( dfx[nX,][,2])
    gtf_file =  paste0( as.vector( dfx[nX,][,3]) ,'/merged.gtf' )

    msg_print( 'Running', nXb )

    # Retrive Data from cufflinks
    R_cuffDB  = readCufflinks( as.character(splr_loc), gtfFile= as.character(gtf_file), genome="mm10", rebuild=TRUE, verbose=TRUE )

    # CummeRbund plots
    R_cuffDB_d = dispersionPlot(genes(R_cuffDB))
    ggsave(R_cuffDB_d ,filename=make_cummeRbund_plot_fp( splr_loc, paste0(nXb,"_cb_dispersion" ) ,".png"))
    R_cuffDB_pBoxRep = csBoxplot(genes(R_cuffDB),replicates=T)
    ggsave(R_cuffDB_pBoxRep ,filename=make_cummeRbund_plot_fp( splr_loc, paste0(nXb,"_cb_pboxrep" ) ,".png"))
    R_cuffDB_pBox = csBoxplot(genes(R_cuffDB))
    ggsave(R_cuffDB_pBox ,filename=make_cummeRbund_plot_fp( splr_loc, paste0(nXb,"_cb_pbox" ) ,".png"))
    R_cuffDB_scatter = csScatterMatrix(genes(R_cuffDB))
    ggsave(R_cuffDB_scatter ,filename=make_cummeRbund_plot_fp( splr_loc, paste0(nXb,"_cb_scatter" ) ,".png"))
    R_cuffDB_genes = fpkmSCVPlot(genes(R_cuffDB))
    ggsave(R_cuffDB_genes ,filename=make_cummeRbund_plot_fp( splr_loc, paste0(nXb,"_cb_fpkm_genes" ) ,".png"))
    R_cuffDB_isoforms = fpkmSCVPlot(isoforms(R_cuffDB))
    ggsave(R_cuffDB_isoforms ,filename=make_cummeRbund_plot_fp( splr_loc, paste0(nXb,"_cb_fpkm_isoforms" ) ,".png"))
    R_cuffDB_densRep = csDensity(genes(R_cuffDB),replicates=T)
    ggsave(R_cuffDB_densRep ,filename=make_cummeRbund_plot_fp( splr_loc, paste0(nXb,"_cb_densnsityrep" ) ,".png"))
    R_cuffDB_mCount = MAplot(genes(R_cuffDB),Cond_1, Cond_2 ,useCount=T)
    ggsave(R_cuffDB_mCount ,filename=make_cummeRbund_plot_fp( splr_loc, paste0(nXb,"_cb_MAplotCount" ) ,".png"))
    R_cuffDB_volcano = csVolcano(genes(R_cuffDB),Cond_1, Cond_2)
    ggsave(R_cuffDB_volcano ,filename=make_cummeRbund_plot_fp( splr_loc, paste0(nXb,"_cb_volcano" ) ,".png"))
    
    # Start spliceR
    R_cuff_spliceR = prepareCuff( R_cuffDB, removeNonChanonicalChr=TRUE )
    spliceRList = spliceR( R_cuff_spliceR, compareTo="preTranscript", filters=filters, expressionCutoff=expressionCutoff, useProgressBar=FALSE )
    #spliceRList = spliceR( R_cuff_spliceR, compareTo="preTranscript", filters=filters, useProgressBar=FALSE )

    # AnnotatePTC , which requires mm10_ucscCDS and save the GTF files
    spliceRList_annoPTC= annotatePTC( spliceRList, mm10_ucscCDS, Mmusculus )
    generateGTF(spliceRList_annoPTC, filters=filters, filePrefix=make_fp( splr_loc, paste0(nXb,"_spliceR" )))
    write.table( as.data.frame( spliceRList_annoPTC$transcript_features ),  make_fp(splr_loc, paste0( nXb , "_spliceRlist_transcriptf_with_annoPTC.csv")),sep=",", quote=FALSE, row.names=FALSE )
    write.table( as.data.frame( spliceRList_annoPTC$exon_features ),  make_fp(splr_loc, paste0( nXb , "_spliceRlist_exonf.csv")),sep=",", quote=FALSE, row.names=FALSE)

    # Run splicR plots and get the data in a data frame
    collect_p = list()
    collect_q = list()
    counter = 1
    for( i in evalv ){
        for ( x in visType ){
            fp_out_pdf = make_spliceR_plot_fp( splr_loc, paste0(nXb,"_spliceR" ), i, x,'.pdf' )
            msg_print('Writing :', fp_out_pdf)
            p = paste0( i,"_",x )
			q = totalNumberOfAS( run_spliceR_plots( spliceRList, evaluate=i, asType=x, out_pdf_fp=fp_out_pdf ) )
			collect_p[[ counter ]] = p
			collect_q[[ counter ]] = q
			counter = counter + 1 
        }
    }
    data_df = data.frame( collect_q )
    colnames( data_df ) = unlist( collect_p )
    data_df  = t( data_df )
    write.table( data_df, make_fp(splr_loc, paste0( nXb , "_spliceRlist_isoform_values.csv")), sep=",", quote=FALSE, row.names=TRUE)
}

print( "done")

