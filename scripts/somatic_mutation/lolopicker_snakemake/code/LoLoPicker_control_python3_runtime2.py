#!/app/easybuild/software/Python/3.7.4-foss-2016b/bin/python -u
import sys, getopt
#import re
import pysam
#import multiprocessing as mp
from multiprocessing import Pool
from datetime import datetime

def main(argv):
    base_quality = 30
    mapping_quality = 30
    thread = 1
    if len(sys.argv) < 3:
        print('usage: LoLoPicker_control.py -l <samplelist> -r <reference> -n <thread> -o <outputpath>')
        sys.exit(1)
    try:
        opts, args = getopt.getopt(argv,"hl:r:n:o:B:m", ["help","samplelist=", "reference=", "thread=", "outputpath=", "basequality=", "mappingquality="])
    except getopt.GetoptError:
        print('usage: LoLoPicker_control.py -l <samplelist> -r <reference> -n <thread> -o <outputpath>')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print('usage: LoLoPicker_control.py -l <samplelist> -r <reference> -n <thread> -o <outputpath>')
            sys.exit()
        elif opt in ("-l", "--samplelist"):
            samplelist = arg
        elif opt in ("-r", "--reference"):
            reference = arg
        elif opt in ("-n", "--thread"):
            thread = arg
        elif opt in ("-o", "--outputpath"):
            outputpath = arg
        elif opt in ("-B", "--basequality"):
            base_quality = arg
        elif opt in ("-m", "--mappingquality"):
            mapping_quality = arg
    return samplelist, reference, thread, outputpath, int(base_quality), int(mapping_quality)

if __name__ == '__main__':
    (samplelist, reference, thread, outputpath, base_quality, mapping_quality) = main(sys.argv[1:])
    ref = pysam.FastaFile(reference)
    ######fixed spelling of "variants" --Anna 11/18/19
    varfile = outputpath + "/raw_somatic_variants.txt"
    tempfile = outputpath + "/control_stats.txt"

    def process_reads(columns, columns_pos, refseq):
        alt_A=0; alt_C=0; alt_G=0; alt_T=0; refcount=0
        ###I think these altReadPos variables can be removed.  --Anna 9/13/2019###
        foundReadName=set(); altReadPos_f_A=[]; altReadPos_r_A=[]; altReadPos_f_G=[]; altReadPos_r_G=[]; altReadPos_f_C=[]; altReadPos_r_C=[]; altReadPos_f_T=[]; altReadPos_r_T=[]
        pileupcolumn = columns;
        if pileupcolumn.pos == columns_pos:
            for pileupread in pileupcolumn.pileups:
                if pileupread.alignment.is_proper_pair and not pileupread.alignment.is_duplicate:
                    if not pileupread.is_del and not pileupread.is_refskip:
                        if (pileupread.alignment.mapping_quality >= mapping_quality) and (pileupread.alignment.query_qualities[pileupread.query_position] >= base_quality):
                            altseq = pileupread.alignment.query_sequence[pileupread.query_position]
                            if pileupread.alignment.query_name not in foundReadName:
                                foundReadName.add(pileupread.alignment.query_name)
                                if pileupread.alignment.query_sequence[pileupread.query_position] == refseq.upper():
                                    refcount += 1
                                else:
                                    if altseq == "A":
                                        alt_A += 1
                                    elif altseq == "G":
                                        alt_G += 1
                                    elif altseq == "C":
                                        alt_C += 1
                                    elif altseq == "T":
                                        alt_T += 1
            return alt_A, alt_C, alt_G, alt_T, refcount

    sample = []
    test_variants = []
    def lookup_controls(test_variants, bampath, sampleID):
        variant_hash = {};
        samfile = pysam.AlignmentFile(bampath)
        print("working on: " + str(sampleID))
        for variant in test_variants:
            #print(variant[0], variant[1], "at", datetime.now(), flush = True) ###Added by Anna 11/18/19 for verification purposes ##flush added 1/8/2020
            reftotal=0;
            chr = variant[0];
            pos = variant[1]-1;
            ref_seq = variant[2];
            alt_base = variant[3]
            for c_columns in samfile.pileup(chr, int(pos), int(pos)+1, truncate=True):
                if int(c_columns.pos) == int(pos):
                    (c_alt_A, c_alt_C, c_alt_G, c_alt_T, c_refcount) = process_reads(c_columns, c_columns.pos, ref_seq)
                    reftotal = reftotal + c_refcount
                    k = (chr+"\t"+str(pos+1)+"\t"+ref_seq+"\t"+alt_base)
                    if alt_base == "A":
                        variant_hash[k] = (str(c_alt_A)+":"+str(c_refcount)+":"+sampleID)
                    if alt_base == "G":
                        variant_hash[k] = (str(c_alt_G)+":"+str(c_refcount)+":"+sampleID)
                    if alt_base == "C":
                        variant_hash[k] = (str(c_alt_C)+":"+str(c_refcount)+":"+sampleID)
                    if alt_base == "T":
                        variant_hash[k] = (str(c_alt_T)+":"+str(c_refcount)+":"+sampleID)
        return variant_hash

    inputlist = open(samplelist)
    for line in ( raw.strip().split() for raw in inputlist):
        sample.append(line)
    inputlist.close()

    vartemp = open(varfile)
    next(vartemp)
    for j in ( raw.strip().split() for raw in vartemp ):
        (chr, pos, ref, alt, t_refcount, t_alt_count, n_refcount, n_alt_count, judge) = j
        if judge == "pass_to_test":
            test_variants.append([str(chr), int(pos), str(ref), str(alt)])
    vartemp.close()

    def ff(sample):
        bampath =  sample[0]; sampleID =  sample[1]
        (control_variants) = lookup_controls(test_variants, bampath, sampleID)
        return control_variants

    p = Pool(int(thread))
    control_variants = p.map(ff, sample)
    merged_variants = {}
    for control_sample in control_variants:
        for k in list(control_sample.keys()):
            if k in merged_variants:
                merged_variants[k] =  merged_variants[k] + "," + control_sample[k]
            else:
                merged_variants[k] = control_sample[k]

    ftemp = open(tempfile, 'w'); vartemp = open(varfile)
    for variants in ( raw.strip().split() for raw in vartemp ):
        var = (variants[0]+"\t"+str(variants[1])+"\t"+variants[2]+"\t"+variants[3])
        if var in merged_variants:
            print('\t'.join(variants)+"\t"+merged_variants[var], file=ftemp)
        else:
            print('\t'.join(variants)+"\t"+"NA", file=ftemp)

    vartemp.close()
    ftemp.close()
    print("Done control panel inspection!", flush = True)
