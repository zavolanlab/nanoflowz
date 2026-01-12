#!/usr/bin/env python3

from __future__ import (absolute_import, division,
                        print_function, unicode_literals)
import warnings
warnings.simplefilter('ignore')

import sys
from argparse import ArgumentParser, RawTextHelpFormatter
import os
import subprocess
import pandas as pd
import numpy as np
import HTSeq


def main():
    """ (1) Extract ONT long reads (Direct RNA sequencing protocol, sequenced from 3' to 5' end) from bam files.
        (2) Create longest possible isoforms based on them as well as annotated transcripts, under assumption that transcript molecules may be partially degraded in 5' to 3' direction AND 3'end cleavage may be not precise (25 nt window).
        (3) Assign isoforms to longest compatible annotated transcripts, under assumption that 5'end of the transcript is well annotated, so isoform 5' end may not exceed that of the annotated isoform.
    """

    __doc__ = "Creating isoforms, annotating them and assigning reads to them"

    parser = ArgumentParser(description=__doc__,
                            formatter_class=RawTextHelpFormatter)

    parser.add_argument("--input_bam_files",
                        dest="input_bam_files",
                        help="Path to the text file, listing the paths to bam files to process",
                        required=True,
                        metavar="FILE",)
    parser.add_argument("--input_gtf_file",
                        dest="input_gtf_file",
                        help="Path to the gtf annotation file",
                        required=True,
                        metavar="FILE",)
    parser.add_argument("--tmp_dir",
                    dest="tmp_dir",
                    help="Path to the temporary dir",
                    required=True,
                    metavar="FILE",)
    parser.add_argument("--out_reads",
                        dest="out_reads",
                        help="path to the output tsv file with reads assigned to isoform ids",
                        required=True,
                        metavar="DIR")
    parser.add_argument("--out_isoforms",
                        dest="out_isoforms",
                        help="path for the output tsv file with assembled isoforms and their annotation",
                        required=True,
                        metavar="FILE")
    parser.add_argument("--ThreePrimeEnd_clustering_distance", 
                        dest="ThreePrimeEnd_clustering_distance",
                        type=int,
                        default=25,
                        help=("integer value for the maximal distance between 3'end sites of reads to be clustered [default: 25]") )
    parser.add_argument("--gtf_skip_rows", 
                        dest="gtf_skip_rows",
                        type=int,
                        default=5,
                        help=("integer value for how many rows to skip when loading gtf file. For GENCODE gtf, usually 5 rows need to be skipped [default: 5]") )    
    try:
        options = parser.parse_args()
    except(Exception):
        parser.print_help()

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    # check that there is '/' sign in the end of dir path
    temp_dir_path = options.tmp_dir
    if len(temp_dir_path)>0:
        if (not temp_dir_path.endswith('/')):
            temp_dir_path = temp_dir_path+'/'
    else:
        print('[FATAL] temp dir is not specified\n')
        return None
    
    almnt_file_paths = pd.read_csv(options.input_bam_files,delimiter="\t",index_col=None,header=None)
    almnt_file_paths = list(almnt_file_paths[0])
    # for further compact storage, assign integer ids to input bam files
    sample_ids = list(range(0,len(almnt_file_paths)))

    print('[INFO] Parsing reads from bam files using HTSeq package\n')
    k = 0
    read_coords = []
    read_introns = []
    for almnt_file_path in almnt_file_paths:
        print('[INFO] processing bam file: '+almnt_file_path+'\n')
        almnt_file = HTSeq.SAM_Reader( almnt_file_path )
        # all .iv elements follow bed format
        j=0
        for almnt in almnt_file:
            pt = almnt.optional_field("pt")
            qs = almnt.optional_field("qs")
            read_coords.append([k,j,almnt.iv.chrom,almnt.iv.start+1,almnt.iv.end,almnt.iv.strand,qs,pt]) # convert start coordinates to 1-based format
            for cigar_elem in almnt.cigar:
                if cigar_elem.type=='N': # splice junction
                    reference_region = cigar_elem.ref_iv
                    read_introns.append([k,j,reference_region.start,reference_region.end+1]) # splice junstion coordinates are defined as last and first 1-based positions of flanking exons
            j=j+1
            if j%10000==0:
                print('[INFO] '+str(j)+' reads done\n')
        if j==0:
            print('[WARNING] No reads in the bam file: '+almnt_file_path+'\n')
        k=k+1
        print('[INFO] Done with bam file: '+almnt_file_path+'\n')

    ###
    # get read information as data frame
    # this is done both for + and - strands
    ###
    read_coords_df = pd.DataFrame(read_coords,columns=['sample_id','within_sample_read_id','chr','start','end','strand','qs','pt'])
    if len(read_coords_df)==0:
        print('[WARNING] No reads extracted. Writing empty output files.\n')
        read_coords_df_final = pd.DataFrame([],columns=['chr', 'start', 'end', 'strand', 'qs', 'pt', 'read_id', 'isoform_id'])
        isoforms_df_final = pd.DataFrame([],columns=['chr', 'strand', 'contig', 'cluster_id', 'isoform_start', 'isoform_end','isoform_id', 'num_of_splice_sites', 'total_transcript_length','annotated_transcript_id'])
        
        read_coords_df_final.to_csv(options.out_reads, sep=str('\t'),header=True,index=None)
        isoforms_df_final.to_csv(options.out_isoforms, sep=str('\t'),header=True,index=None)        
        return None
    
    read_coords_df['read_id'] = read_coords_df['sample_id'].astype('str')+';'+read_coords_df['within_sample_read_id'].astype('str')
    
    ###
    # get split-read information as data frame
    # this is done both for + and - strands
    ###
    
    read_introns_df = pd.DataFrame(read_introns,columns=['sample_id','within_sample_read_id','start','end'])
    read_introns_df['read_id'] = read_introns_df['sample_id'].astype('str')+';'+read_introns_df['within_sample_read_id'].astype('str')
    
    # print important info
    print('[INFO] Summary stats about reads from bam files:\n')
    k=0
    for almnt_file_path in almnt_file_paths:
        info = almnt_file_path+'\nnumber of reads: '+str(len(read_coords_df.loc[read_coords_df['sample_id']==k]))
        info = info+'\n'
        info = info+'among them, number of split reads: '+str(len(read_introns_df.loc[read_introns_df['sample_id']==k]['within_sample_read_id'].unique()))
        print(info)
        k=k+1 
    
    print('\n\
    #################\n\
    [INFO] We then merged reads from samples\n')
    
    read_coords_df = read_coords_df.drop(['sample_id','within_sample_read_id'],axis=1)
    read_introns_df = read_introns_df.drop(['sample_id','within_sample_read_id'],axis=1)
    read_introns_df = pd.merge(read_introns_df,read_coords_df[['read_id','strand','chr']],how='left',on='read_id')
    print('[INFO] Total number of reads: '+str(len(read_coords_df))+'\n')

    print('\n\
    #################\n\
    [INFO] Parsing GTF file\n')
    
    # get annotation file and extract annotated transcript isoforms in the same manner as observed ones
    # here, 1-based coordinates are expected
    gtf = pd.read_csv(options.input_gtf_file,delimiter="\t",index_col=None,header=None,skiprows=options.gtf_skip_rows)

    if len(gtf)==0:
        print('[FATAL] GTF file is empty\n')
        return None
    
    # extract transcript elements
    transcripts = gtf.loc[gtf[2]=='transcript'].reset_index(drop=True)
    transcripts['transcript_id'] = transcripts[8].str.split('transcript_id "',expand=True)[1].str.split('";',expand=True)[0]
    transcript_ends_bed = pd.concat([transcripts.loc[transcripts[6]=='+'][[0,4,'transcript_id',6]].reset_index(drop=True),
                               transcripts.loc[transcripts[6]=='-'][[0,3,'transcript_id',6]].reset_index(drop=True).rename(columns={3:4})]).reset_index(drop=True)
    transcript_ends_bed.columns = ['chr','end','read_id','strand']

    # extract exon elements
    exons = gtf.loc[gtf[2]=='exon'].reset_index(drop=True)
    exons['transcript_id'] = exons[8].str.split('transcript_id "',expand=True)[1].str.split('";',expand=True)[0]

    # obtain annotated introns on + strand
    exons_plus = exons.loc[exons[6]=='+'].sort_values([0,'transcript_id',3,4],ascending=[True,True,True,True]).reset_index(drop=True)
    cur_transcript = ''
    a = []
    for elem in exons_plus[[0,3,4,6,'transcript_id']].values:
        if elem[4]==cur_transcript:
            a.append([elem[0],pred_end,elem[1],elem[3],cur_transcript])
        cur_transcript = elem[4]
        pred_end = elem[2]
    introns_plus = pd.DataFrame(a,
                                columns=['chr','start','end','strand','read_id']).sort_values(
        ['chr','strand','read_id','end','start'],ascending=[True,True,True,False,False]).reset_index(drop=True)
    introns_plus['coords'] = introns_plus['end'].astype('str')+','+introns_plus['start'].astype('str')
    print('[INFO] '+str(len(introns_plus[['chr','strand','coords']].drop_duplicates()))+' introns in '+str(len(introns_plus['read_id'].unique()))+' transcripts were extracted from gtf on + strand\n')
    introns_plus = pd.concat([introns_plus[['read_id','chr','strand']],introns_plus.groupby('read_id')['coords'].transform(lambda x: ','.join(x))],axis=1).drop_duplicates().reset_index(drop=True)
    introns_plus = pd.merge(introns_plus,
                            transcripts[['transcript_id',3,4]].rename(columns={'transcript_id':'read_id',3:'start',4:'end'}),
                            how='left',on='read_id')
    
    # obtain annotated introns on - strand
    exons_minus = exons.loc[exons[6]=='-'].sort_values([0,'transcript_id',3,4],ascending=[True,True,True,True]).reset_index(drop=True)
    cur_transcript = ''
    a = []
    for elem in exons_minus[[0,3,4,6,'transcript_id']].values:
        if elem[4]==cur_transcript:
            a.append([elem[0],pred_end,elem[1],elem[3],cur_transcript])
        cur_transcript = elem[4]
        pred_end = elem[2]
    introns_minus = pd.DataFrame(a,
                                columns=['chr','start','end','strand','read_id']).sort_values(
        ['chr','strand','read_id','start','end'],ascending=[True,True,True,True,True]).reset_index(drop=True)
    introns_minus['coords'] = introns_minus['end'].astype('str')+','+introns_minus['start'].astype('str')
    print('[INFO] '+str(len(introns_minus[['chr','strand','coords']].drop_duplicates()))+' introns in '+str(len(introns_minus['read_id'].unique()))+' transcripts were extracted from gtf on - strand\n')
    introns_minus = pd.concat([introns_minus[['read_id','chr','strand']],introns_minus.groupby('read_id')['coords'].transform(lambda x: ','.join(x))],axis=1).drop_duplicates().reset_index(drop=True)
    introns_minus = pd.merge(introns_minus,
                            transcripts[['transcript_id',3,4]].rename(columns={'transcript_id':'read_id',3:'start',4:'end'}),
                            how='left',on='read_id')
    

    ###
    # first, cluster the 3'end coordinates of reads
    # this is done both for + and - strands
    # bedtools is used for that
    ###
    
    read_ends_bed = pd.concat([read_coords_df.loc[read_coords_df['strand']=='+'][['chr','end','read_id','strand']].reset_index(drop=True),
                               read_coords_df.loc[read_coords_df['strand']=='-'][['chr','start','read_id','strand']].reset_index(drop=True).rename(columns={'start':'end'})]).reset_index(drop=True)
    # add annotated isoforms
    read_ends_bed = pd.concat([read_ends_bed,transcript_ends_bed]).reset_index(drop=True)
    
    read_ends_bed['start'] = read_ends_bed['end']-1
    read_ends_bed['score'] = 0
    read_ends_bed[['chr','start','end','read_id','score','strand']].to_csv(temp_dir_path+'read_ends.bed', sep=str('\t'),header=False,index=None)
    
    print("""[INFO] sorting 3'ends of reads with bedtools sort\n""")
    command = 'bedtools sort -i '+temp_dir_path+'read_ends.bed'+' > '+temp_dir_path+'read_ends.sorted.bed'
    out = subprocess.check_output(command, shell=True)
    
    print("""[INFO] clustering 3'ends of reads with bedtools cluster\n""")
    command = 'bedtools cluster -d '+str(options.ThreePrimeEnd_clustering_distance)+' -s -i '+temp_dir_path+'read_ends.sorted.bed > '+temp_dir_path+'read_ends.clustered.bed'
    out = subprocess.check_output(command, shell=True)
    
    read_ends_clustered = pd.read_csv(temp_dir_path+'read_ends.clustered.bed',delimiter="\t",index_col=None,header=None)
    
    if len(read_ends_bed)>0 and len(read_ends_clustered)>0 and len(read_ends_bed)==len(read_ends_clustered):
        print("""[INFO] clustering 3'ends of reads was successful\n""")
    else:
        print("""[FATAL] clustering 3'ends of reads was NOT successful\n""")
        return None
    
    ###
    # Define splice isoforms, based on split reads and annotated transcripts with introns
    # isoforms are defined by the cluser_id of 3'end, combination of splice sites (contig), and most upstream start site
    ###
    
    ###
    ## process positive (+) strand
    ###
    
    print('\n\
    #################\n\
    [INFO] process positive (+) strand\n\
    #################\n')
    
    # sort splits in reads from closest to 3' to farthest from 3', because reads are sequenced in 3' -> 5' direction, 
    # so 5' end of the read is not informative, it may not necessarily reflect true transcript end, rather just sequencing abruption, or also RNA degradation
    tmp_plus = read_introns_df.loc[read_introns_df['strand']=='+'].sort_values(['read_id','end','start'],ascending=[True,False,False]).reset_index(drop=True)
    tmp_plus['coords'] = tmp_plus['end'].astype('str')+','+tmp_plus['start'].astype('str')
    # now, each row corresponds to one intron in the read
    # get the sequence of splice sites within reads (those are expected to be 1-based coordinates of exon ends)
    tmp_plus = pd.concat([tmp_plus[['read_id','chr','strand']],tmp_plus.groupby('read_id')['coords'].transform(lambda x: ','.join(x))],axis=1).drop_duplicates().reset_index(drop=True)
    # now, each row corresponds to one read
    tmp_plus = pd.merge(tmp_plus,read_coords_df[['read_id','start','end']],how='left',on='read_id')
    tmp_plus = pd.merge(tmp_plus,read_ends_clustered[[3,6]].rename(columns = {3:'read_id',6:'cluster_id'}),how='left',on='read_id') # add cluster ids of 3'ends of reads
    
    if 'cluster_id' in list(introns_plus.columns):
        introns_plus = introns_plus.drop(['cluster_id'],axis=1)
    introns_plus = pd.merge(introns_plus,read_ends_clustered[[3,6]].rename(columns = {3:'read_id',6:'cluster_id'}),how='left',on='read_id') # add cluster ids of 3'ends of reads
    
    # append annotated introns on + strand, but only from the transcripts that have the same cluster ids as actual split reads 
    tmp_plus = pd.concat([tmp_plus,
                          introns_plus.loc[introns_plus['cluster_id'].isin(list(tmp_plus['cluster_id'].unique()))]]).reset_index(drop=True) 
    
    
    tmp_plus = tmp_plus.sort_values(['chr','strand','cluster_id','coords']).reset_index(drop=True) # coords are sorted lexicographically
    
    # create isoforms, taking advantage of lexicographic sorting
    
    cur_cluster,cur_contig,contigs = -1,'',[] 
    # here, "contig" is referred to the longest combination of splice sites that fully contains one or more preceeding (in terms of sorting) combinations
    for elem in tmp_plus[['coords','cluster_id']].values:
        if (not elem[0].startswith(cur_contig)) or (elem[1]!=cur_cluster and cur_cluster!=(-1)):
            contigs.append([cur_contig,cur_cluster])
        cur_contig,cur_cluster = elem[0],elem[1]
    contigs.append([cur_contig,cur_cluster])
    contigs_df = pd.DataFrame(contigs,columns = ['contig','cluster_id'])
    isoforms_df = pd.merge(contigs_df,
            tmp_plus[['coords','cluster_id','chr','strand','start','end']].rename(columns={'coords':'contig'}),
            how='left',
            on=['contig','cluster_id']).groupby(['chr','strand','contig','cluster_id']).agg({'start':np.min,'end':np.max}).reset_index()
    isoforms_df[['start','end']] = isoforms_df[['start','end']].astype('int')
    isoforms_df['isoform_id'] = (isoforms_df.index).astype('str')+'_+'
    isoforms_df = isoforms_df.rename(columns={'start':'isoform_start','end':'isoform_end'})
    isoforms_df['num_of_splice_sites'] = isoforms_df.apply(lambda x:len(x['contig'].split(',')),1)
    
    def get_transcript_length(x):
        total_len = 0
        prev_pos = x['isoform_end']
        i=0
        for pos in x['contig'].split(','):
            if i==0:
                total_len = total_len+(prev_pos-int(pos)+1)
            elif i%2==0:
                total_len = total_len+(prev_pos-int(pos)+1)
            prev_pos = int(pos)
            i=i+1
        total_len = total_len+(prev_pos-x['isoform_start']+1)
        return total_len
    isoforms_df['total_transcript_length'] = isoforms_df.apply(lambda x:get_transcript_length(x),1)
    
    tmp_plus = pd.merge(tmp_plus,isoforms_df,how='left',on=['chr','strand','cluster_id']) # (1) match reads to isoforms by cluser_id of 3'end
    tmp_plus['splice_compatible'] = tmp_plus.apply(lambda x:x['contig'].startswith(x['coords']),1) 
    tmp_plus_sel = tmp_plus.loc[tmp_plus['splice_compatible']].reset_index(drop=True) # (2) match reads to isoforms by splice contig compatibility
    
    # (3) finally, assign to longest compatible isoform; if several isoforms meet the criteria, randomly assign to one of them
    read_to_isoform_assignment = pd.merge(tmp_plus_sel,
             tmp_plus_sel.groupby(['read_id']).agg({'total_transcript_length':max}).reset_index(),
            how='inner',
             on=['read_id','total_transcript_length'])[['read_id','isoform_id']].drop_duplicates(['read_id']).reset_index(drop=True)

    
    isoform_annotated = pd.merge(read_to_isoform_assignment.loc[read_to_isoform_assignment['read_id'].str.startswith('ENSMUS')],
                                tmp_plus_sel[['read_id','start','end','coords']],
                                how='left',
                                on = ['read_id'])

    if len(isoform_annotated)>0:
        isoform_annotated = pd.merge(isoform_annotated,isoforms_df[['isoform_id','isoform_start','isoform_end']],how='left',on='isoform_id')
        # require that the start of annotated isoform is equal or longer relative to isoform start 
        isoform_annotated = isoform_annotated.loc[isoform_annotated['start']<=isoform_annotated['isoform_start']].reset_index(drop = True)
        
        # to resolve possible duplicates, select the longest compatible annotated transcript
        def get_transcript_length(x):
            total_len = 0
            prev_pos = x['end']
            i=0
            for pos in x['coords'].split(','):
                if i==0:
                    total_len = total_len+(prev_pos-int(pos)+1)
                elif i%2==0:
                    total_len = total_len+(prev_pos-int(pos)+1)
                prev_pos = int(pos)
                i=i+1
            total_len = total_len+(prev_pos-x['start']+1)
            return total_len
        isoform_annotated['total_transcript_length'] = isoform_annotated.apply(lambda x:get_transcript_length(x),1)
        isoform_annotated = pd.merge(isoform_annotated,
                                    isoform_annotated.groupby(['isoform_id']).agg({'total_transcript_length':max}).reset_index(),
                                    how='inner',on=['isoform_id','total_transcript_length']).drop_duplicates('isoform_id').reset_index(drop=True)
        
        # add annotated transcripts that correspond to the isoform
        isoforms_df = pd.merge(isoforms_df,
                    isoform_annotated[['isoform_id','read_id']].rename(columns={'read_id':'annotated_transcript_id'}),
                    how = 'left',
                    on='isoform_id')
        isoforms_df['annotated_transcript_id'] = isoforms_df['annotated_transcript_id'].fillna('')
    else:
        isoforms_df['annotated_transcript_id'] = ''
        
    # remove annotated transcripts from data frame with isoform assignment
    read_to_isoform_assignment = read_to_isoform_assignment.loc[~read_to_isoform_assignment['read_id'].str.startswith('ENSMUS')].reset_index(drop=True)
    
    # fraction of isoforms without assigned transcript
    print('[INFO] Fraction of spliced isoforms without assigned transcript: '+str(len(isoforms_df.loc[isoforms_df['annotated_transcript_id']==''])/len(isoforms_df)))
    
    # check that there are no transcripts to which multimple isoforms are assigned
    isoforms_df['t']=1
    gr = isoforms_df.groupby('annotated_transcript_id').agg({'t':sum}).reset_index()
    print('[INFO] max number of transcripts assigned to an isoform: '+str(gr.loc[gr['annotated_transcript_id']!='']['t'].max()))
    if gr.loc[gr['annotated_transcript_id']!='']['t'].max()>1:
        print('[FATAL] Max number of transcripts assigned to an isoform should be 1\n')
        return None

    isoforms_df = isoforms_df.drop(['t'],axis=1)
    
    # take reads on + strand that were not assigned yet -> try to assign them to compatible splice isoforms 
    # (under assumption that those reads are due to sequencing abruption / RNA degradation in 5'->3' direction)
    
    tmp_plus_rest = read_coords_df.loc[(read_coords_df['strand']=='+')&(~read_coords_df['read_id'].isin(list(read_to_isoform_assignment['read_id'])))].reset_index(drop=True)
    tmp_plus_rest = pd.merge(tmp_plus_rest,read_ends_clustered[[3,6]].rename(columns = {3:'read_id',6:'cluster_id'}),how='left',on='read_id') # add cluster ids of 3'ends of reads
    tmp_plus_rest = pd.merge(tmp_plus_rest,isoforms_df,how='left',on=['chr','strand','cluster_id']) # (1) match reads to existing splice isoforms by cluser_id of 3'end
    tmp_plus_rest = tmp_plus_rest.loc[~tmp_plus_rest['contig'].isna()].reset_index(drop=True)
    
    tmp_plus_rest['terminal_exon_start_coord'] = tmp_plus_rest.apply(lambda x:int(x['contig'].split(',')[0]),1)
    tmp_plus_rest = tmp_plus_rest.loc[tmp_plus_rest['start']>=tmp_plus_rest['terminal_exon_start_coord']].reset_index(drop=True) # require that reads do not contradict the splice contig
    tmp_plus_rest['total_transcript_length'] = tmp_plus_rest['total_transcript_length'].astype('int')
    read_to_isoform_assignment_chunk = pd.merge(tmp_plus_rest,
                        tmp_plus_rest.groupby(['read_id']).agg({'total_transcript_length':np.max}).reset_index(),
                        how='inner',
                         on=['read_id','total_transcript_length'])[['read_id','isoform_id']].drop_duplicates(['read_id']).reset_index(drop=True)
    read_to_isoform_assignment = pd.concat([read_to_isoform_assignment,read_to_isoform_assignment_chunk]).reset_index(drop=True) # append to assigned reads
    
    # take reads on + strand that were still not assigned yet
    tmp_plus_rest = read_coords_df.loc[(read_coords_df['strand']=='+')&(
        ~read_coords_df['read_id'].isin(list(read_to_isoform_assignment['read_id'])))].reset_index(drop=True)[['chr','start','end','strand','read_id']]
    tmp_plus_rest = pd.merge(tmp_plus_rest,read_ends_clustered[[3,6]].rename(columns = {3:'read_id',6:'cluster_id'}),how='left',on='read_id') # add cluster ids of 3'ends of reads
    # take single-exon annotated transcripts
    single_exon_transcripts = transcripts.loc[(transcripts[6]=='+')&(transcripts['transcript_id'].isin(list(introns_plus['read_id'].unique())))][[0,3,4,6,'transcript_id']]
    single_exon_transcripts.columns = ['chr','start','end','strand','read_id']
    single_exon_transcripts = pd.merge(single_exon_transcripts,read_ends_clustered[[3,6]].rename(columns = {3:'read_id',6:'cluster_id'}),how='left',on='read_id') # add cluster ids of 3'ends of reads
    # leave only those annotated single-exon transcripts whose cluster ids (3'ends) are among those that reads have
    single_exon_transcripts = single_exon_transcripts.loc[single_exon_transcripts['cluster_id'].isin(list(tmp_plus_rest['cluster_id'].unique()))].reset_index(drop=True)
    # make a union of reads with single-exon transcripts
    tmp_plus_rest = pd.concat([tmp_plus_rest,single_exon_transcripts]).reset_index(drop=True)
    
    # create isoforms
    isoforms_df_chunk = tmp_plus_rest.groupby(['chr','strand','cluster_id']).agg({'start':np.min,'end':np.max}).reset_index().rename(columns={'start':'isoform_start','end':'isoform_end'})
    isoforms_df_chunk[['isoform_start','isoform_end']] = isoforms_df_chunk[['isoform_start','isoform_end']].astype('int')
    isoforms_df_chunk['contig'] = ''
    isoforms_df_chunk['num_of_splice_sites'] = 0
    isoforms_df_chunk['total_transcript_length'] = (isoforms_df_chunk['isoform_end']-isoforms_df_chunk['isoform_start']+1)
    isoforms_df_chunk['isoform_id'] = (max(isoforms_df.index)+1+isoforms_df_chunk.index).astype('str')+'_+'
    
    isoform_annotated = pd.merge(tmp_plus_rest.loc[tmp_plus_rest['read_id'].str.startswith('ENSMUS')],
            isoforms_df_chunk[['cluster_id','isoform_start','isoform_end','isoform_id']],
            how='left',on='cluster_id')
    isoform_annotated = isoform_annotated.loc[isoform_annotated['start']<=isoform_annotated['isoform_start']].reset_index(drop=True)
    isoform_annotated['total_transcript_length'] = isoform_annotated['end']-isoform_annotated['start']+1
    isoform_annotated['total_transcript_length'] = isoform_annotated['total_transcript_length'].astype('int')
    isoform_annotated = pd.merge(isoform_annotated,
                                isoform_annotated.groupby(['isoform_id']).agg({'total_transcript_length':np.max}).reset_index(),
                                how='inner',on=['isoform_id','total_transcript_length']).drop_duplicates('isoform_id').reset_index(drop=True)
    
    # add annotated transcripts that correspond to the isoform
    isoforms_df_chunk = pd.merge(isoforms_df_chunk,
                isoform_annotated[['isoform_id','read_id']].rename(columns={'read_id':'annotated_transcript_id'}),
                how = 'left',
                on='isoform_id')
    isoforms_df_chunk['annotated_transcript_id'] = isoforms_df_chunk['annotated_transcript_id'].fillna('')
    
    print('[INFO] Fraction of single-exon isoforms without an assigned transcript: '+str(len(isoforms_df_chunk.loc[isoforms_df_chunk['annotated_transcript_id']==''])/len(isoforms_df_chunk)))
    
    isoforms_df = pd.concat([isoforms_df,isoforms_df_chunk[list(isoforms_df.columns)]]).reset_index(drop=True) # append to existing ones
    
    read_to_isoform_assignment_chunk = pd.merge(tmp_plus_rest.loc[~tmp_plus_rest['read_id'].str.startswith('ENSMUS')],
                                                isoforms_df,
                                                how='left',
                                                on=['chr','strand','cluster_id'])[['read_id','isoform_id']].drop_duplicates(['read_id']).reset_index(drop=True)
    read_to_isoform_assignment_plus = pd.concat([read_to_isoform_assignment,read_to_isoform_assignment_chunk]).reset_index(drop=True) # append to assigned reads
    # here, for final assingment, separate plus and minus strands
    
    print('[INFO] Number of reads on + strand without assigned isoform: '+str(len(read_coords_df.loc[(read_coords_df['strand']=='+')&(
        ~read_coords_df['read_id'].isin(list(read_to_isoform_assignment_plus['read_id'])))].reset_index(drop=True))))
    
    # only keep isoforms that were assigned to at least one read
    isoforms_df_plus = isoforms_df.loc[isoforms_df['isoform_id'].isin(list(read_to_isoform_assignment_plus['isoform_id'].unique()))].reset_index(drop=True)
    
    print('[INFO] Number of: \n\
    (1) reads on + strand, \n\
    (2) reads on + strand with assigned isoform, \n\
    (3) assigned isoform ids for reads on + strand, \n\
    (4) retrieved isoforms on + strand: \n'+str([len(read_coords_df.loc[(read_coords_df['strand']=='+')]),
        len(read_to_isoform_assignment_plus['read_id'].unique()),
        len(read_to_isoform_assignment_plus['isoform_id'].unique()),
        len(isoforms_df_plus)]))

    if len(read_coords_df.loc[(read_coords_df['strand']=='+')])!=len(read_to_isoform_assignment_plus['read_id'].unique()):
        print('[FATAL] Number (1) is not equal to number (2)\n')
        return None
    if len(read_to_isoform_assignment_plus['isoform_id'].unique())!=len(isoforms_df_plus):
        print('[FATAL] Number (3) is not equal to number (4)\n')
        return None
    
    print('\n\
    #################\n\
    [INFO] process negative (-) strand\n\
    #################\n')
    
    # sort splits in reads from closest to 3' to farthest from 3', because reads are sequenced in 3' -> 5' direction, 
    # so 5' end of the read is not informative, it may not necessarily reflect true transcript end, rather just sequencing abruption, or also RNA degradation
    tmp_minus = read_introns_df.loc[read_introns_df['strand']=='-'].sort_values(['read_id','start','end'],ascending=[True,True,True]).reset_index(drop=True)
    tmp_minus['coords'] = tmp_minus['start'].astype('str')+','+tmp_minus['end'].astype('str')
    # now, each row corresponds to one intron in the read
    # get the sequence of splice sites within reads (those are expected to be 1-based coordinates of exon ends)
    tmp_minus = pd.concat([tmp_minus[['read_id','chr','strand']],tmp_minus.groupby('read_id')['coords'].transform(lambda x: ','.join(x))],axis=1).drop_duplicates().reset_index(drop=True)
    # now, each row corresponds to one read
    tmp_minus = pd.merge(tmp_minus,read_coords_df[['read_id','start','end']],how='left',on='read_id')
    tmp_minus = pd.merge(tmp_minus,read_ends_clustered[[3,6]].rename(columns = {3:'read_id',6:'cluster_id'}),how='left',on='read_id') # add cluster ids of 3'ends of reads
    
    if 'cluster_id' in list(introns_minus.columns):
        introns_minus = introns_minus.drop(['cluster_id'],axis=1)
    introns_minus = pd.merge(introns_minus,read_ends_clustered[[3,6]].rename(columns = {3:'read_id',6:'cluster_id'}),how='left',on='read_id') # add cluster ids of 3'ends of reads
    
    # append annotated introns on - strand, but only from the transcripts that have the same cluster ids as actual split reads 
    tmp_minus = pd.concat([tmp_minus,
                          introns_minus.loc[introns_minus['cluster_id'].isin(list(tmp_minus['cluster_id'].unique()))]]).reset_index(drop=True) 
    
    tmp_minus = tmp_minus.sort_values(['chr','strand','cluster_id','coords']).reset_index(drop=True) # coords are sorted lexicographically
    
    # create isoforms, taking advantage of lexicographic sorting
    
    cur_cluster,cur_contig,contigs = -1,'',[] 
    # here, "contig" is referred to the longest combination of splice sites that fully contains one or more preceeding (in terms of sorting) combinations
    for elem in tmp_minus[['coords','cluster_id']].values:
        if (not elem[0].startswith(cur_contig)) or (elem[1]!=cur_cluster and cur_cluster!=(-1)):
            contigs.append([cur_contig,cur_cluster])
        cur_contig,cur_cluster = elem[0],elem[1]
    contigs.append([cur_contig,cur_cluster])
    contigs_df = pd.DataFrame(contigs,columns = ['contig','cluster_id'])
    
    isoforms_df = pd.merge(contigs_df,
            tmp_minus[['coords','cluster_id','chr','strand','start','end']].rename(columns={'coords':'contig'}),
            how='left',
            on=['contig','cluster_id']).groupby(['chr','strand','contig','cluster_id']).agg({'start':np.min,'end':np.max}).reset_index()
    isoforms_df[['start','end']] = isoforms_df[['start','end']].astype('int')
    isoforms_df['isoform_id'] = (isoforms_df.index).astype('str')+'_-'
    isoforms_df = isoforms_df.rename(columns={'start':'isoform_start','end':'isoform_end'})
    isoforms_df['num_of_splice_sites'] = isoforms_df.apply(lambda x:len(x['contig'].split(',')),1)
    
    def get_transcript_length(x):
        total_len = 0
        prev_pos = x['isoform_start']
        i=0
        for pos in x['contig'].split(','):
            if i==0 or i%2==0:
                total_len = total_len+(int(pos)-prev_pos+1)
            prev_pos = int(pos)
            i=i+1
        total_len = total_len+(x['isoform_end']-prev_pos+1)
        return total_len
    isoforms_df['total_transcript_length'] = isoforms_df.apply(lambda x:get_transcript_length(x),1)
    
    tmp_minus = pd.merge(tmp_minus,isoforms_df,how='left',on=['chr','strand','cluster_id']) # (1) match reads to isoforms by cluser_id of 3'end
    tmp_minus['splice_compatible'] = tmp_minus.apply(lambda x:x['contig'].startswith(x['coords']),1)
    
    tmp_minus_sel = tmp_minus.loc[tmp_minus['splice_compatible']].reset_index(drop=True) # (2) match reads to isoforms by splice contig compatibility
    
    # (3) finally, assign to longest compatible isoform; if several isoforms meet the criteria, randomly assign to one of them
    tmp_minus_sel['total_transcript_length'] = tmp_minus_sel['total_transcript_length'].astype('int')
    read_to_isoform_assignment = pd.merge(tmp_minus_sel,
             tmp_minus_sel.groupby(['read_id']).agg({'total_transcript_length':np.max}).reset_index(),
            how='inner',
             on=['read_id','total_transcript_length'])[['read_id','isoform_id']].drop_duplicates(['read_id']).reset_index(drop=True)
    
    isoform_annotated = pd.merge(read_to_isoform_assignment.loc[read_to_isoform_assignment['read_id'].str.startswith('ENSMUS')],
                                tmp_minus_sel[['read_id','start','end','coords']],
                                how='left',
                                on = ['read_id'])

    if len(isoform_annotated)>0:    
        isoform_annotated = pd.merge(isoform_annotated,isoforms_df[['isoform_id','isoform_start','isoform_end']],how='left',on='isoform_id')
        # require that the end of annotated isoform is equal or longer relative to isoform end 
        isoform_annotated = isoform_annotated.loc[isoform_annotated['end']>=isoform_annotated['isoform_end']].reset_index(drop = True)
        
        # to resolve possible duplicates, select the longest compatible annotated transcript
        def get_transcript_length(x):
            total_len = 0
            prev_pos = x['start']
            i=0
            for pos in x['coords'].split(','):
                if i==0 or i%2==0:
                    total_len = total_len+(int(pos)-prev_pos+1)
                prev_pos = int(pos)
                i=i+1
            total_len = total_len+(x['end']-prev_pos+1)
            return total_len
        isoform_annotated['total_transcript_length'] = isoform_annotated.apply(lambda x:get_transcript_length(x),1).astype('int')
        isoform_annotated = pd.merge(isoform_annotated,
                                    isoform_annotated.groupby(['isoform_id']).agg({'total_transcript_length':np.max}).reset_index(),
                                    how='inner',on=['isoform_id','total_transcript_length']).drop_duplicates('isoform_id').reset_index(drop=True)
        
        # add annotated transcripts that correspond to the isoform
        isoforms_df = pd.merge(isoforms_df,
                    isoform_annotated[['isoform_id','read_id']].rename(columns={'read_id':'annotated_transcript_id'}),
                    how = 'left',
                    on='isoform_id')
        isoforms_df['annotated_transcript_id'] = isoforms_df['annotated_transcript_id'].fillna('')
    else:
        isoforms_df['annotated_transcript_id'] = ''
        
    # remove annotated transcripts from data frame with isoform assignment
    read_to_isoform_assignment = read_to_isoform_assignment.loc[~read_to_isoform_assignment['read_id'].str.startswith('ENSMUS')].reset_index(drop=True)
    
    # fraction of isoforms without assigned transcript
    print('[INFO] Fraction of spliced isoforms without assigned transcript: '+str(len(isoforms_df.loc[isoforms_df['annotated_transcript_id']==''])/len(isoforms_df)))
    
    # check that there are no transcripts to which multimple isoforms are assigned
    isoforms_df['t']=1
    gr = isoforms_df.groupby('annotated_transcript_id').agg({'t':sum}).reset_index()
    print('[INFO] max number of transcripts assigned to an isoform: '+str(gr.loc[gr['annotated_transcript_id']!='']['t'].max()))
    if gr.loc[gr['annotated_transcript_id']!='']['t'].max()>1:
        print('[FATAL] Max number of transcripts assigned to an isoform should be 1\n')
        return None
    
    isoforms_df = isoforms_df.drop(['t'],axis=1)
    
    # take reads on - strand that were not assigned yet -> try to assign them to compatible splice isoforms 
    # (under assumption that those reads are due to sequencing abruption / RNA degradation in 5'->3' direction)
    
    tmp_minus_rest = read_coords_df.loc[(read_coords_df['strand']=='-')&(~read_coords_df['read_id'].isin(list(read_to_isoform_assignment['read_id'])))].reset_index(drop=True)
    tmp_minus_rest = pd.merge(tmp_minus_rest,read_ends_clustered[[3,6]].rename(columns = {3:'read_id',6:'cluster_id'}),how='left',on='read_id') # add cluster ids of 3'ends of reads
    tmp_minus_rest = pd.merge(tmp_minus_rest,isoforms_df,how='left',on=['chr','strand','cluster_id']) # (1) match reads to existing splice isoforms by cluser_id of 3'end
    tmp_minus_rest = tmp_minus_rest.loc[~tmp_minus_rest['contig'].isna()].reset_index(drop=True)
    
    tmp_minus_rest['terminal_exon_end_coord'] = tmp_minus_rest.apply(lambda x:int(x['contig'].split(',')[0]),1)
    tmp_minus_rest = tmp_minus_rest.loc[tmp_minus_rest['end']<=tmp_minus_rest['terminal_exon_end_coord']].reset_index(drop=True) # require that reads do not contradict the splice contig
    tmp_minus_rest['total_transcript_length'] = tmp_minus_rest['total_transcript_length'].astype('int')
    read_to_isoform_assignment_chunk = pd.merge(tmp_minus_rest,
                        tmp_minus_rest.groupby(['read_id']).agg({'total_transcript_length':np.max}).reset_index(),
                        how='inner',
                         on=['read_id','total_transcript_length'])[['read_id','isoform_id']].drop_duplicates(['read_id']).reset_index(drop=True) # assign to longest compatible splice isoform
    read_to_isoform_assignment = pd.concat([read_to_isoform_assignment,read_to_isoform_assignment_chunk]).reset_index(drop=True) # append to assigned reads
    
    # take reads on - strand that were still not assigned yet
    tmp_minus_rest = read_coords_df.loc[(read_coords_df['strand']=='-')&(
        ~read_coords_df['read_id'].isin(list(read_to_isoform_assignment['read_id'])))].reset_index(drop=True)[['chr','start','end','strand','read_id']]
    tmp_minus_rest = pd.merge(tmp_minus_rest,read_ends_clustered[[3,6]].rename(columns = {3:'read_id',6:'cluster_id'}),how='left',on='read_id') # add cluster ids of 3'ends of reads
    # take single-exon annotated transcripts
    single_exon_transcripts = transcripts.loc[(transcripts[6]=='-')&(transcripts['transcript_id'].isin(list(introns_minus['read_id'].unique())))][[0,3,4,6,'transcript_id']]
    single_exon_transcripts.columns = ['chr','start','end','strand','read_id']
    single_exon_transcripts = pd.merge(single_exon_transcripts,read_ends_clustered[[3,6]].rename(columns = {3:'read_id',6:'cluster_id'}),how='left',on='read_id') # add cluster ids of 3'ends of reads
    # leave only those annotated single-exon transcripts whose cluster ids (3'ends) are among those that reads have
    single_exon_transcripts = single_exon_transcripts.loc[single_exon_transcripts['cluster_id'].isin(list(tmp_minus_rest['cluster_id'].unique()))].reset_index(drop=True)
    # make a union of reads with single-exon transcripts
    tmp_minus_rest = pd.concat([tmp_minus_rest,single_exon_transcripts]).reset_index(drop=True)
    
    # create isoforms
    isoforms_df_chunk = tmp_minus_rest.groupby(['chr','strand','cluster_id']).agg({'start':np.min,'end':np.max}).reset_index().rename(columns={'start':'isoform_start','end':'isoform_end'})
    isoforms_df_chunk['contig'] = ''
    isoforms_df_chunk['num_of_splice_sites'] = 0
    isoforms_df_chunk['total_transcript_length'] = (isoforms_df_chunk['isoform_end']-isoforms_df_chunk['isoform_start']+1).astype('int')
    isoforms_df_chunk['isoform_id'] = (np.max(isoforms_df.index)+1+isoforms_df_chunk.index).astype('str')+'_-'
    
    isoform_annotated = pd.merge(tmp_minus_rest.loc[tmp_minus_rest['read_id'].str.startswith('ENSMUS')],
            isoforms_df_chunk[['cluster_id','isoform_start','isoform_end','isoform_id']],
            how='left',on='cluster_id')
    isoform_annotated = isoform_annotated.loc[isoform_annotated['end']>=isoform_annotated['isoform_end']].reset_index(drop=True)
    isoform_annotated['total_transcript_length'] = (isoform_annotated['end']-isoform_annotated['start']+1).astype('int')
    isoform_annotated = pd.merge(isoform_annotated,
                                isoform_annotated.groupby(['isoform_id']).agg({'total_transcript_length':np.max}).reset_index(),
                                how='inner',on=['isoform_id','total_transcript_length']).drop_duplicates('isoform_id').reset_index(drop=True)
    
    # add annotated transcripts that correspond to the isoform
    isoforms_df_chunk = pd.merge(isoforms_df_chunk,
                isoform_annotated[['isoform_id','read_id']].rename(columns={'read_id':'annotated_transcript_id'}),
                how = 'left',
                on='isoform_id')
    isoforms_df_chunk['annotated_transcript_id'] = isoforms_df_chunk['annotated_transcript_id'].fillna('')
    
    print('[INFO] Fraction of single-exon isoforms without an assigned transcript: '+str(len(isoforms_df_chunk.loc[isoforms_df_chunk['annotated_transcript_id']==''])/len(isoforms_df_chunk)))
    
    isoforms_df = pd.concat([isoforms_df,isoforms_df_chunk[list(isoforms_df.columns)]]).reset_index(drop=True) # append to existing ones
    
    read_to_isoform_assignment_chunk = pd.merge(tmp_minus_rest.loc[~tmp_minus_rest['read_id'].str.startswith('ENSMUS')],
                                                isoforms_df,
                                                how='left',
                                                on=['chr','strand','cluster_id'])[['read_id','isoform_id']].drop_duplicates(['read_id']).reset_index(drop=True)
    read_to_isoform_assignment_minus = pd.concat([read_to_isoform_assignment,read_to_isoform_assignment_chunk]).reset_index(drop=True) # append to assigned reads
    # here, for final assingment, separate plus and minus strands
    
    print('[INFO] Number of reads on - strand without assigned isoform: '+str(len(read_coords_df.loc[(read_coords_df['strand']=='-')&(
        ~read_coords_df['read_id'].isin(list(read_to_isoform_assignment_minus['read_id'])))].reset_index(drop=True))))
    
    # only keep isoforms that were assigned to at least one read
    isoforms_df_minus = isoforms_df.loc[isoforms_df['isoform_id'].isin(list(read_to_isoform_assignment_minus['isoform_id'].unique()))].reset_index(drop=True)
    
    print('[INFO] Number of: \n\
    (1) reads on - strand, \n\
    (2) reads on - strand with assigned isoform, \n\
    (3) assigned isoform ids for reads on - strand, \n\
    (4) retrieved isoforms on - strand: \n'+str([len(read_coords_df.loc[(read_coords_df['strand']=='-')]),
        len(read_to_isoform_assignment_minus['read_id'].unique()),
        len(read_to_isoform_assignment_minus['isoform_id'].unique()),
        len(isoforms_df_minus)]))

    if len(read_coords_df.loc[(read_coords_df['strand']=='-')])!=len(read_to_isoform_assignment_minus['read_id'].unique()):
        print('[FATAL] Number (1) is not equal to number (2)\n')
        return None
    if len(read_to_isoform_assignment_minus['isoform_id'].unique())!=len(isoforms_df_minus):
        print('[FATAL] Number (3) is not equal to number (4)\n')
        return None    
    
    read_coords_df_final = pd.merge(read_coords_df,
                pd.concat([read_to_isoform_assignment_plus,read_to_isoform_assignment_minus]).reset_index(drop=True),
                 how='left',
                 on='read_id')
    isoforms_df_final = pd.concat([isoforms_df_plus,isoforms_df_minus]).reset_index(drop=True)    

    print('[INFO] Writing output files\n')
    read_coords_df_final.to_csv(options.out_reads, sep=str('\t'),header=True,index=None)
    isoforms_df_final.to_csv(options.out_isoforms, sep=str('\t'),header=True,index=None)
    
    
if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        sys.stderr.write("User interrupt!")
        sys.exit(1)