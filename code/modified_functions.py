import numpy as np
import pandas as pd
import argparse
from tqdm import tqdm
import gzip
from gtfparse import read_gtf

# def write_annotation(gtf, out_file, gene_type = 'protein_coding'):
    
#     if gene_type != 'all':
#         print('Working on ' + gene_type + ' genes')
        
#         try:
#             gene_list = gtf.loc[gtf.gene_type == gene_type].gene.unique()
#         except:
#             raise Exception(
#                 'Gene type not found in GTF file. Please make sure that gene_type annotation exists, or use "--gene_type all" instead.'
#             )
        
#     else:
#         print('Working on all genes')
#         gene_list = gtf.gene.unique()
        
#     fh = gzip.open(out_file + '.tab.gz', 'wt')
    
#     fh.write('\t'.join(['', 'intron', 'event', 'gene', 'transcript', 'exon_number']) + '\n')
    
#     gene_counts = 1
#     total_genes = len(gene_list)
    
#     print('Working in a total of ' + str(total_genes))
    
#     for gene in tqdm(gene_list, leave=True, position=0):
#         rows = process_gene_annotation(gtf, gene)
#         print(rows)
        
#         fh.write(rows)
            
#         gene_counts += 1
        
#     print('Processed ' + str(gene_counts-1)+'/'+str(total_genes))
    
#     print('Getting constitutive introns')
#     get_const_introns(fh, gtf)
    
#     fh.close()
#     print('Finished writing annotation')

def write_annotation(gtf, out_file, gene_type='protein_coding'):
    
    if gene_type != 'all':
        print('Working on ' + gene_type + ' genes')
        try:
            gene_list = gtf.loc[gtf.gene_type == gene_type].gene.unique()
        except:
            raise Exception('Gene type not found in GTF file.')
    else:
        print('Working on all genes')
        gene_list = gtf.gene.unique()
    
    with gzip.open(out_file + '.tab.gz', 'wt') as fh:  # context manager guarantees flush+close
        fh.write('\t'.join(['', 'intron', 'event', 'gene', 'transcript', 'exon_number']) + '\n')

        total_genes = len(gene_list)
        print('Working on a total of ' + str(total_genes))

        for i, gene in enumerate(tqdm(gene_list, leave=True, position=0)):
            rows = process_gene_annotation(gtf, gene)
            fh.write(rows)
            if i % 100 == 0:
                fh.flush()
            
    print('Getting constitutive introns') 
    with gzip.open(out_file + '.tab.gz', 'at') as fh:
        get_const_introns(fh, gtf)
    
    print('Finished writing annotation')

def get_cassette_exons(gtf_gene, donor_dir, acceptor_dir):

    ase_dir = {}

    for donor in donor_dir.keys():
        if len(donor_dir[donor].keys()) > 1:                    # Check for donors with more than one acceptor
            acceptors = sorted(donor_dir[donor].keys())         # Sort them, so that we go in order for acceptors

            for n in range(1, len(acceptors)):                     # start with the 2nd acceptor; we assume no exon inside 1st intron

                i = acceptors[n]                                   # parse through acceptors
                intron_list = donor_dir[donor][i]                  # introns that end in the nth acceptor
                for intron in intron_list:                      # parsing through all the introns

                    for j in acceptor_dir[i].keys():            # for each donor of the nth acceptor

                        if j == donor:
                            skipped_isoforms = [x for x in donor_dir[j][i] if x in acceptor_dir[i][j]] 

                        if j != donor:                          # check one that is different

                            cassette_isoforms = []
                            for acceptor_isoform in acceptor_dir[i][j]:
                                split_intron = acceptor_isoform.split('.')

                                upstream_intron = split_intron[0] + '.' + str(int(split_intron[1]) - 1)

                                for acceptor_final in donor_dir[donor].keys():
                                    if upstream_intron in donor_dir[donor][acceptor_final]:
                                        I1 = donor + ':' + acceptor_final
                                        I2 = j + ':' + i
                                        SE = donor + ':' + i

                                        evento = (I1, I2, SE)

                                        # cassette_isoforms.append(upstream_intron + ':' + acceptor_isoform)

                                        nmd = upstream_intron.split('.')[0]
                                        
                                        try:

                                            nmd_type = gtf_gene.loc[gtf_gene.transcript == nmd].transcript_type.unique()[0]
                                            
                                        except:
                                            nmd_type = 'notype'
                                    
                                        event_key = '|'.join(evento)
                                        # Get exon # for particular isoform
                                        exon_number = gtf_gene.loc[
                                            (gtf_gene.transcript == nmd) & 
                                            (gtf_gene.start == (int(acceptor_final) + 1)) & 
                                            (gtf_gene.end == (int(j) - 1)),
                                            'exon_number'
                                        ]
                                        transcript_id = upstream_intron.split('.')[0]  # same as nmd, already computed
                                        entry = (transcript_id, nmd_type, int(exon_number))
                                        if event_key not in ase_dir:
                                            ase_dir[event_key] = [entry]
                                        else:
                                            ase_dir[event_key].append(entry)
                                        
    return ase_dir

def process_gene_annotation(gtf, gene):
    
    gtf_gene = gtf.loc[gtf.gene == gene]
    
    donor_dir, acceptor_dir = get_dirs(gtf_gene, gene)
    
    ase_dir = get_cassette_exons(gtf_gene, donor_dir, acceptor_dir)
    
    chrom = list(gtf_gene.chrom)[0]
    strand = list(gtf_gene.strand)[0]
    
    event_sorted = sorted(ase_dir.keys())
    
    gene_events = ''

    notype_counts = 1
    pc_counts = 1
    nmd_counts = 1
    other_counts = 1
    
    for event in event_sorted:
        
        for transcript, transcript_type, exon_number in ase_dir[event]:
            if 'protein_coding' in transcript_type:
                event_name = gene + '_ProteinCoding_' + str(pc_counts)
                pc_counts += 1
            elif 'nonsense_mediated_decay' in transcript_type:
                event_name = gene + '_NMD_' + str(nmd_counts)
                nmd_counts += 1
            ### modify the script here if you're interested in other tags
            ### elif 'your_tag' in ase_dir[event]:
            ###    event_name = gene + '_YOURTAG_' + str(yourtag_counts)
            elif 'notype' in transcript_type:
                event_name = gene + '_UNKNOWNTYPE_' + str(notype_counts)
                notype_counts += 1
            else:
                event_name = gene + '_other_' + str(other_counts)
                other_counts += 1
                
            i1_start, i1_end = event.split('|')[0].split(':')
            intron_loc = chrom + ':' + i1_start + '-' + i1_end + ':' + strand
            i1_line = '\t'.join([event_name + '_I1', intron_loc, event_name, gene, transcript, str(exon_number)]) + '\n'
            gene_events += i1_line
            
            
            i2_start, i2_end = event.split('|')[1].split(':')
            intron_loc = chrom + ':' + i2_start + '-' + i2_end + ':' + strand
            i2_line = '\t'.join([event_name + '_I2', intron_loc, event_name, gene, transcript, str(exon_number)]) + '\n'
            gene_events += i2_line
            
            se_start, se_end = event.split('|')[2].split(':')
            intron_loc = chrom + ':' + se_start + '-' + se_end + ':' + strand
            se_line = '\t'.join([event_name + '_SE', intron_loc, event_name, gene, transcript, str(exon_number)]) + '\n'
            gene_events += se_line        

    return gene_events
       
  
    
def process_gtf(gtf_file, exclude, gene_name, no_trim_id, gene_type_tag, transcript_type_tag):

    print('Processing GTF file...')
    
    exclude_chromosomes = exclude.split(',')
    
    gtf = read_gtf(gtf_file)
    gtf = gtf.to_pandas()

    try:

        try:
            cols = ['seqname', 'start', 'end', 'strand', 'transcript_id', gene_type_tag, gene_name, transcript_type_tag, 'exon_number']
            gtf = gtf.loc[pd.Series([x not in exclude_chromosomes for x in gtf.seqname]) & (gtf.feature == 'exon'), cols]
            gtf.columns = ['chrom', 'start', 'end', 'strand', 'transcript', 'gene_type', 'gene', 'transcript_type', 'exon_number']
        except:
            cols = ['seqname', 'start', 'end', 'strand', 'transcript_id', gene_name, 'exon_number']
            gtf = gtf.loc[pd.Series([x not in exclude_chromosomes for x in gtf.seqname]) & (gtf.feature == 'exon'), cols]
            gtf.columns = ['chrom', 'start', 'end', 'strand', 'transcript', 'gene', 'exon_number']
            
    except:
        raise Exception('Isufficient information to create annotation. transcript_id is needed to find cassette exons.')

    if not no_trim_id:
        gtf['gene'] = gtf.gene.str.replace(r'\.\d+$', '', regex=True)

    gtf.transcript = [x.split('.')[0] for x in gtf.transcript]

    gtf = gtf.loc[gtf.gene != '']
        
    return gtf

#######


def get_dirs(gtf_gene, gene):
    donor_dir = {}
    acceptor_dir = {}

    old_transcript = ''
    for idx, row in gtf_gene.sort_values(['transcript', 'start', 'end']).iterrows():
        transcript = row.transcript

        if old_transcript != transcript:

            counter = 1
            old_transcript = transcript

        else:

            transcript_e = transcript + '.' + str(counter)
            acceptor = str(int(row.start)-1)

            if not acceptor in acceptor_dir.keys():
                acceptor_dir.update({acceptor:{donor:[transcript_e]}})

            elif not donor in acceptor_dir[acceptor].keys():
                acceptor_dir[acceptor].update({donor:[transcript_e]})

            else:
                acceptor_dir[acceptor][donor].append(transcript_e)


            if acceptor in donor_dir[donor].keys():
                donor_dir[donor][acceptor].append(transcript_e)

            else:
                donor_dir[donor].update({acceptor:[transcript_e]})

            counter += 1

        donor = str(int(row.end)+1)

        if not donor in donor_dir.keys():
            #print('repetido')
            donor_dir.update({donor:{}})

    return donor_dir, acceptor_dir

    
def get_const_introns(fh, gtf):
    
    try:
        pc_gtf = gtf.loc[(gtf.transcript_type == 'protein_coding') & (gtf.gene_type == 'protein_coding')]
    except:
        pc_gtf = gtf
    
    for i in tqdm(pc_gtf.groupby('gene'), leave=True, position=0):
        t_gtf = i[1]
        transcripts = t_gtf.transcript.unique()

        donor_dir, acceptor_dir = get_dirs(t_gtf, transcripts[0])

        introns_list = []
        for donor in donor_dir:

            if len(donor_dir[donor].keys()) > 0:

                acceptor = list(donor_dir[donor].keys())[0]

                intron = str(t_gtf.chrom.unique()[0])
                intron += ':' + str(donor)
                intron += ':' + str(acceptor)
                intron += ':' + str(t_gtf.strand.unique()[0])

                introns_list.append(intron)

        for t in transcripts:

            t_gtf = i[1].loc[i[1].transcript == t]

            introns = []

            donor_dir, acceptor_dir = get_dirs(t_gtf, transcripts[0])

            for donor in donor_dir:

                if len(donor_dir[donor].keys()) > 0:

                    acceptor = list(donor_dir[donor].keys())[0]

                    intron = str(t_gtf.chrom.unique()[0])
                    intron += ':' + str(donor)
                    intron += ':' + str(acceptor)
                    intron += ':' + str(t_gtf.strand.unique()[0])

                    introns.append(intron)

            introns_list = [x for x in introns_list if x in introns]

        intron_counts = 1
        if len(introns_list) > 0:
            for intron in introns_list:
                chrom, start, end, strand = intron.split(':')

                intron_number = i[0] + '_' + str(intron_counts)
                intron_name = intron_number + '_CI'
                intron_loc = chrom + ':' + start + '-' + end + ':' + strand
                
                linea = '\t'.join([intron_name, intron_loc, intron_number, i[0], '', '']) + '\n'
                
                fh.write(linea)

                intron_counts += 1

    