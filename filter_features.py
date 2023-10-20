import argparse
import subprocess
import os
import random
import pysam
import sys
import numpy as np

parser = argparse.ArgumentParser(
                    prog = 'Filter features',
                    description = 'Script for removing reads that map within given annotation features in case of multi-mapping reads, if one of the mapping position is within teh feature, remove read completely with all other alignments',
                    epilog = '...')

#parser.add_argument('mode', type=str, help='Transcriptomic or genomic mode')
parser.add_argument('-f','--features', type=str, required=True, dest='features', nargs='+', help='Feature that will be filtered. mRNA, tRNA or rRNA')
parser.add_argument('-gf', '--gtf', type=str, required=True, dest='gtf', help='Path to annotation file in gtf format')
parser.add_argument('-b', '--bam', type=str, required=True, dest='bam', help='Path to alignments file in sorted bam format')
parser.add_argument('-o', '--out_bam', type=str, default='filtered_features.bam', dest='out_bam', help='Path to output bam file in which filtered alignments will be stored')
parser.add_argument('-r', '--reject', action='store_true', dest='reject', help='If used reads mapped to filtered feaature will be removed')
parser.add_argument('-w', '--write', action='store_true', help='If used output BAM file will be saved. Default - do not write BAM file, only statistics on stdout')

args = parser.parse_args()
r = str(int(random.random() * 1000000))
print(r)

for feature in args.features:
    print(feature)
    out_name = "{}_filtered_{}.gtf.bed".format(r,feature)
    tmp_file = '{}_tmp_reads.txt'.format(r)
    if feature == 'mRNA':
        cmd = 'grep -P "^\\S+\\s+\\S+\\s+transcript" {} | grep  "gbkey \\\"Gene\\\"" | grep "gene_biotype \\\"protein_coding\\\""  | cut -f 1,4,5 > {}'.format(args.gtf,out_name)
    else:
        cmd = 'grep -P "^\\S+\\s+\\S+\\s+transcript" {} | grep  "gbkey \\"{}\\"" | cut -f 1,4,5 > {}'.format(args.gtf,feature,out_name)
    os.system(cmd)
    print('Filtering features done')
    cmd = 'bedtools intersect -bed -abam {} -b {} | cut -f 4'.format(args.bam, out_name)
    filtered_reads = subprocess.getoutput(cmd).split('\n')
    print('Filtering reads done')
    cmd = 'samtools view -F 260 {} | cut -f 1'.format(args.bam)
    reads = subprocess.getoutput(cmd).split('\n')
    print('Subtracting primary alignments done')
    intersection = np.intersect1d(np.array(reads), np.array(filtered_reads))
    inter_list = intersection.tolist()
    print("partially done")
    print(tmp_file)
    out_tmp = open(tmp_file,'w')
    s = '\n'.join(inter_list)
    out_tmp.write(s)
    out_tmp.close()
    print('Intercextion between all and filtered reads done')
    if args.write:
        if args.out_bam == 'filtered_features.bam':
            args.out_bam = '{}_filtered_features.bam'.format(r)
        #samfile = pysam.AlignmentFile(args.bam, 'rb')
        #out = pysam.AlignmentFile(args.out_bam, 'wb', header=samfile.header)
        if args.reject:
            cmd = 'samtools view -h {} | grep -v -f {} > {}'.format(args.bam, tmp_file, args.out_bam)
            os.system(cmd)
            '''for read in samfile:
                if read.query_name not in inter_list:
                    out.write(read)'''
        else:
            cmd = 'samtools view -h {} | grep -f {} > {}'.format(args.bam, tmp_file, args.out_bam)
            os.system(cmd)
            '''for read in samfile:
                if read.query_name in inter_list:
                    out.write(read)'''
        #samfile.close()
        #out.close()
        print('New BAM done')
    overlapping_reads = len(intersection)
    all_reads = len(reads)
    others_reads = all_reads - overlapping_reads

    print("Features used for filtering: {}".format(feature))
    print("Number of input reads (primary aligned): {}".format(all_reads))
    if (not args.reject):
        print("Number of reatined reads: {}".format(overlapping_reads))
        print("Number of removed reads: {}".format(others_reads))
        print("Percent of retained reads: {}".format(round((overlapping_reads*100)/all_reads,2)))
    else:
        print("Number of reatined reads: {}".format(others_reads))
        print("Number of removed reads: {}".format(overlapping_reads))
        print("Percent of retained reads: {}%").format(round((others_reads*100)/all_reads,2))