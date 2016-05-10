#!/usr/bin/python
import argparse
import re
import subprocess
import sys
import os
import time
import string
import random
__author__ = 'Robin van der Weide, Bsc. (robinweide@gmail.com)'
start_time = time.time()
parser = argparse.ArgumentParser(
    description='This is REPPER: the fast rDNArep-counter.')
parser.add_argument('-bam', '--inputbam', help='BAM file of fastq(s) to REPPER.fa, sorted and filtered (-F4 -q1).', required=False)
parser.add_argument('-fqf', '--inputfastq_forward', help='Fastq-files to be aligned.', required=False)
parser.add_argument('-fqr', '--inputfastq_reverse', help='Fastq-files to be aligned.', required=False)
parser.add_argument('-fq', '--inputfastq', help='Single-end Fastq-file to be aligned.', required=False)
#parser.add_argument('-p', '--ploidy', help='Ploidy of sample(chr21).', required=True)
parser.add_argument('-wb', '--windowedGCBed', help='Windowed bed-file of REPPER.fa, incl gc-percentage per window.', required=False)
parser.add_argument('-w', '--window', help='Selected window.', required=False)
#parser.add_argument('-s', '--slow', help='Slow or fast settings (Y/N).', required=False)
parser.add_argument('-d', '--output', help='Full path of dir (a new folder will be created in this).', required=True)
args = vars(parser.parse_args())

#create job-ID
ID = ''.join(random.choice(string.ascii_uppercase + string.digits) for _ in range(10))

#check is output-dir exists: yes, make sub-dir
if os.path.isdir(args['output']) is not True:
    sys.stderr.write(str("[error:repper]\tOutput-dir ") + str(args['output-dir']) + "does not exist: exiting." +  str(".\n"))
    sys.exit()
else:
    os.makedirs(args['output'] + ID)
    Odir = str(args['output'] + ID + "/")

#function for running shell commmands
def cmdline(command):
    process = subprocess.Popen(
        args=command,
        stdout=subprocess.PIPE,
        shell=True
    )
    return process.communicate()[0]

#print welcome message
sys.stderr.write("\nThis is REPPER: the fast(ish) rDNArep-pipeline." + str("\n\n"))
sys.stderr.write("[info:repper]\tJob-ID is: " + str(ID) + str(".\n"))

#parse name of fastq (either PE or SE)
if args['inputfastq_reverse'] is not None:
    base, tail = str(args['inputfastq_forward']).split('.')
    FASTQ = base
else:
    base, tail = str(args['inputfastq']).split('.')
    FASTQ = base


#change dir to repper-dir
os.chdir(os.path.dirname(os.path.abspath(__file__)))
os.chdir("..")

''' depricated
# choose the fast or slow reference
REF=''
slowARG = ''
totalSeqLenHG = float(0)
if args['slow'] is str("Y"):
    REF = str("REPPER")
    slowARG = str("yes")
    totalSeqLenHG = float(2372470)
else:
    REF = str("REPPER350bp")
    slowARG = str("no")
    totalSeqLenHG = float(1072695)
'''
#set reference values
REF = str("REPPER_V3")
totalSeqLenHG = float(31117029)

#make dict of sequence sizes
sizeFile = open("./data/REPPER_V3.sizes",'r')
totalSeqLenRDNA = float(7663)
sizedict = {}
for row in sizeFile:
     Seq, Size = re.split(r'\t+', row.rstrip())
     SEQ = Seq.upper()
     sizedict[SEQ.rstrip()] = float(Size.rstrip())


#IF BAM IS NOT GIVEN, Do the mapping of a fastq onto REPPER.fa
BAM=''
if args['inputbam'] is None:
    #paired-end fastqs?
    if args['inputfastq_reverse'] is not None:
        sys.stderr.write("[info:mapping]\tMapping " + str(FASTQ) + "." + str("\n"))
        sys.stderr.write("[info:mapping]\tThis might take a long time (WGS 10X with 3 threads: ~4 hours)!" + str("\n"))
        cmd1 = """bwa mem -v 1 -t 3 """
        cmda = os.getcwd()
        cmd2 = """.fa """
        cmd3 = os.path.abspath(args['inputfastq_forward'])
        cmd3b = """ """
        cmd3c = os.path.abspath(args['inputfastq_reverse'])
        cmd4 = """ | samtools view -b -f 2 - | samtools sort - """
        cmda = os.getcwd()
        cmd5 = """_vs_"""
        cmd6 = os.path.splitext(os.path.basename(FASTQ))[0]
        cmd = str(cmd1 + cmda + "/data/" + REF + cmd2 + cmd3 + cmd3b + cmd3c + cmd4 + Odir + REF + cmd5 + cmd6)
        subprocess.call(cmd, shell=True)
        BAM = str(Odir + REF + cmd5 + cmd6 + ".bam")
        sys.stderr.write("[info:mapping]\tMapping of " + str(FASTQ) + " done." + str("\n"))
    #single-end fastq?
    elif args['inputfastq'] is not None:
        sys.stderr.write("[info:mapping]\tMapping " + str(FASTQ) + "." + str("\n"))
        sys.stderr.write("[info:mapping]\tThis might take a long time (WGS 10X with 3 threads: ~4 hours)!" + str("\n"))
        cmd1 = """bwa mem -v 1 -t 3 """
        cmda = os.getcwd()
        cmd2 = """.fa """
        cmd3 = os.path.abspath(args['inputfastq'])
        cmd4 = """ | samtools view -b -F 4 - | samtools sort - """
        cmda = os.getcwd()
        cmd5 = """_vs_"""
        cmd6 = os.path.splitext(os.path.basename(FASTQ))[0]
        cmd = str(cmd1 + cmda + "/data/" + REF + cmd2 + cmd3 + cmd4 + Odir + REF + cmd5 + cmd6)
        subprocess.call(cmd, shell=True)
        BAM = str(Odir + REF + cmd5 + cmd6 + ".bam")
        sys.stderr.write("[info:mapping]\tMapping of " + str(FASTQ) + " done." + str("\n"))
    else:
        sys.stderr.write("[err:mapping]\tNo BAM-file or fastq-file given." + str("\n"))
        sys.stderr.write("[err:mapping]\tGiving up." + str("\n"))
        sys.exit()
else:
    BAM = args['inputbam']

#Make and Load bed-file, containing the desired windows and the gc content
#[contigname     start    stop   ContigWindow#   GC]
windowedGCBedName = ''
if args['windowedGCBed'] is None:
    if args['window'] is not None:
        sys.stderr.write("[info:window]\tNo windows-bed given, but window-parameter does exist." + str("\n"))
        sys.stderr.write("[info:window]\tMaking new windows-bed..." + str("\n"))
        cmd1 = """samtools faidx """
        cmd1b = """.fa; bedtools makewindows -b ."""
        cmd1c = """.bed -i winnum -w """
        cmd2 = """ | bedtools nuc -fi """
        cmd2a = """.fa  -bed stdin | grep -v ^# | awk '{print $1"\\t"$2"\\t"$3"\\t"$4"\\t"$6}' | tee tmp1 | awk '{printf("%.2f\\n", $5)}' > tmp2 ; awk '{print $1"\\t"$2"\\t"$3"\\t"$4}' tmp1 | paste - tmp2 > """
        cmd2c = """_"""
        cmd3 = """.bed ; rm tmp*"""""
        cmd = str(cmd1 + os.getcwd() + "/data/" + REF + cmd1b + "/data/" + REF + cmd1c + args['window'] + cmd2 + os.getcwd() + "/data/" + REF + cmd2a + Odir + REF + cmd2c + args['window'] + cmd3)
        subprocess.call(cmd, shell=True)
        windowedGCBedName = str(Odir + REF + "_" + args['window'] + ".bed")
        windowedGCBed = open(str(Odir + REF + "_" + args['window'] + ".bed"), 'r')
        sys.stderr.write("[info:window]\tDone." + str("\n"))
    else:
        args['window'] = float(5000)
        sys.stderr.write("[info:window]\tNo windows-bed given." + str("\n"))
        windowedGCBed = open(str(str(".") + "/data/" + REF + str("_5000.bed")), 'r')
        windowedGCBedName = str(str(".") + "/data/" + REF + str("_5000.bed"))
        sys.stderr.write("[info:window]\tDefault windowed bed-file (5000bp) loaded." + str("\n"))

else:
    windowedGCBed = open(args['windowedGCBed'], 'r')
    windowedGCBedName = str(args['windowedGCBed'])


#Get per-windows coverage
sys.stderr.write("[info:depth]\tGetting average coverage for each window..." + str("\n"))
bedTEMP = str(Odir + "tmp_" + ID + ".bed")
lines_in_bed = sum(1 for line in windowedGCBed)
refcount = float(0)
WGClist = []

#per-contig coverage
BamBedCMDa = """bedtools genomecov -bga -ibam """
BamBedCMD = str(str(BamBedCMDa) + BAM + str(" > " + Odir + ID + str(".bed")))
cmdline(BamBedCMD)
tmpfile = str(Odir + str(ID) + str(".tmp"))
covbed = str(Odir + ID + str(".bed"))

#per window, get mean coverage
for row in open(windowedGCBedName):
    refcount += float(1)
    chr,start,stop,window,gc = re.split(r'\t+', row.rstrip())
    fo = open(bedTEMP, 'w')
    fo.write(row)
    fo.close()
    wc_2 = """bedtools intersect -a """
    wc_2a = """ -b """
    wc_3 = """ | awk '{ sum += $4 } END {print sum / (0.0001+NR) }' > """
    wc_3a = """ ; awk '{print}' """
    wc_4 = """ | awk '{print $1"_"$4"\\t"$5"\\t"$3-$2}' | paste - """
    tab = """\\t"""
    wc_cmd = str(wc_2 + covbed + wc_2a + bedTEMP + wc_3 + tmpfile + wc_3a + bedTEMP + wc_4 + tmpfile)
    p = cmdline(wc_cmd).strip().lower()
    WGClist.append(p)
    #progress bar 
    perc = float(100) * float(refcount)/float(lines_in_bed)
    if float(perc) < float(100):
        sys.stderr.write("\t\r%    d%%    " % perc)
        sys.stderr.flush()
    else:
        sys.stderr.write("\t\r%    d%%    " % int(100))
        sys.stderr.flush()
sys.stderr.write("\n[info:depth]\tDone." + str("\n"))
cmdline(str("rm ") + covbed)
cmdline(str("rm ") + tmpfile)


#compute avg cov per gc
sys.stderr.write("[info:rmodel]\tComputing regression-curve for LOESS-smoothing." + str("\n"))
GCdict = {}
GCmdict = {}
#for each GC, make a list with coverages
for i in WGClist:
    window,gc,size,cov = re.split(r'\t+', i.rstrip())
    if gc in GCdict:
        GCdict[gc].append(float(cov))
    else:
        GCdict[gc] = [float(cov)]

#for each GC-coverage-list, calculate mean coverage and store in GCmdict
for g,cl in GCdict.items():
    GCmdict[g] = sum(cl) / float(len(cl))
sys.stderr.write("[info:rmodel]\tDone." + str("\n"))
sys.stderr.write("[info:LOESS]\tStarting GC-normalisation by LOESS-smoothing." + str("\n"))

#get average coverage
covList =[]
loesList = []
for i in WGClist:
    window,gc,size,cov = re.split(r'\t+', i.rstrip())
    covList.append(float(cov))
averageCoverage = sum(covList) / float(len(covList))

#Perform the LOESS correction
ro = open(str(Odir + ID + str('_rawOutput')), 'a')
for i in WGClist:
    window,gc,size,cov = re.split(r'\t+', i.rstrip())
    if float(cov) > float(0):
        LOESS = float(cov) - (float(GCmdict[gc])-averageCoverage)
        if LOESS <= float(0):
            outline = str(str(window) + str("\t") + str(gc) + str("\t") + str(cov) + str("\t") + str(float(0)))
            ro.write(outline + str("\n"))
            loesList.append(outline)
        else:
            outline = str(str(window) + str("\t") + str(gc) + str("\t") + str(cov) + str("\t") + str(LOESS))
            ro.write(outline + str("\n"))
            loesList.append(outline)
    else:
        LOESS = float(0)
        outline = str(str(window) + str("\t") + str(gc) + str("\t") + str(cov) + str("\t") + str(LOESS))
        ro.write(outline + str("\n"))
        loesList.append(outline)
ro.close()
sys.stderr.write("[info:LOESS]\tSaved " + str(Odir + ID + str('_rawOutput') + str("\n")))
sys.stderr.write("[info:LOESS]\tDone." + str("\n"))



'''
calculate average of windows for both control and rDNA (subunits)

'''
sys.stderr.write("[info:copy]\tComputing copy-numbers." + str("\n"))

covHGlist = []
covRDNAlist = []
cov28list = []
cov18list = []
cov58list = []
Wcov28list = []
Wcov18list = []
Wcov58list = []
WcovRDNAlist = []
WcovHGlist = []


for row in loesList:
    window,gc,cov,LOESS = re.split(r'\t+', row.rstrip())
    if window.startswith('control'):
        #seqq = "_".join(window.split("_")[:-1])
        #seq = seqq.upper()
        #size = sizedict[seq]
        #windowNo = "_".join(window.split("_")[-1:])
        #maxSize = float(args['window'])*float(windowNo)
        #windowSize = float(0)
        #if float(maxSize) > float(size):
        #    a = float(size)-float(maxSize)
        #    windowSize = float(args['window']) - abs(a)
        #else:
        #    windowSize = float(args['window'])
        #weightedLOESS = float(LOESS)*(windowSize/totalSeqLenHG)
        #WcovHGlist.append(float(weightedLOESS))                        # corrected for length and GC
        covHGlist.append(float(LOESS))                                 # corrected for GC only
    elif window.startswith('rdna'):
        #seqq = "_".join(window.split("_")[:-1])
        #seq = seqq.upper()
        #size = sizedict[seq]
        #windowNo = "_".join(window.split("_")[-1:])
        #maxSize = float(args['window'])*float(windowNo)
        #windowSize = float(0)
        #if float(maxSize) > float(size):
        #    a = float(size)-float(maxSize)
        #    windowSize = float(args['window']) - abs(a)
        #else:
        #    windowSize = float(args['window'])
        #weightedLOESS = float(LOESS)*(windowSize/totalSeqLenRDNA)  
        #WcovRDNAlist.append(float(weightedLOESS))                      # corrected for length and GC
        covRDNAlist.append(float(LOESS))                               # corrected for GC only
        if window.startswith('rdna18s'):
        #    seqq = window[:-3]
        #    seq = seqq.upper()
        #    weightedLOESS = float(LOESS)*(sizedict[seq]/totalSeqLenRDNA)
        #    Wcov18list.append(float(weightedLOESS))                     # corrected for length and GC
            cov18list.append(float(LOESS))                              # corrected for GC only
        elif window.startswith('rdna5_8s'):
        #    seqq = window[:-3]
        #    seq = seqq.upper()
        #    weightedLOESS = float(LOESS)*(sizedict[seq]/totalSeqLenRDNA)
        #    Wcov58list.append(float(weightedLOESS))                     # corrected for length and GC
            cov58list.append(float(LOESS))                              # corrected for GC only
        elif window.startswith('rdna28s'):
        #    seqq = window[:-2]
        #    seq = seqq.upper()
        #    weightedLOESS = float(LOESS)*(sizedict[seq]/totalSeqLenRDNA)
        #    Wcov28list.append(float(weightedLOESS))                     # corrected for length and GC
                cov28list.append(float(LOESS))                              # corrected for GC only
#WmeanCovHG = float(sum(WcovHGlist) / float(len(WcovHGlist)))
#WmeanCovRDNA =  float(sum(WcovRDNAlist) / float(len(WcovRDNAlist)))
#WmeanCov18 =  float(sum(Wcov18list) / float(len(Wcov18list)))
#WmeanCov28 =  float(sum(Wcov28list) / float(len(Wcov28list)))
#WmeanCov58 =  float(sum(Wcov58list) / float(len(Wcov58list)))
meanCovHG = float(sum(covHGlist) / float(len(covHGlist)))
meanCovRDNA =  float(sum(covRDNAlist) / float(len(covRDNAlist)))
meanCov18 =  float(sum(cov18list) / float(len(cov18list)))
meanCov28 =  float(sum(cov28list) / float(len(cov28list)))
meanCov58 =  float(sum(cov58list) / float(len(cov58list)))

copyR = float(meanCovRDNA/meanCovHG)
copy18 = float(meanCov18/meanCovHG)
copy28 = float(meanCov28/meanCovHG)
copy58 = float(meanCov58/meanCovHG)

#WcopyR = float(WmeanCovRDNA/WmeanCovHG)
#Wcopy18 = float(WmeanCov18/WmeanCovHG)
#Wcopy28 = float(WmeanCov28/WmeanCovHG)
#Wcopy58 = float(WmeanCov58/WmeanCovHG)

#sys.stderr.write("[result:copy]\tMean coverage rDNA: " + str(meanCovRDNA)+ str("\n"))
#sys.stderr.write("[result:copy]\tMean coverage 18s: " + str(meanCov18)+ str("\n"))
#sys.stderr.write("[result:copy]\tMean coverage 5.8s: " + str(meanCov58)+ str("\n"))
#sys.stderr.write("[result:copy]\tMean coverage 28s: " + str(meanCov28)+ str("\n"))
#sys.stderr.write("[result:copy]\tMean coverage hg: " + str(meanCovHG)+ str("\n"))
#sys.stderr.write("[result:copy]\tCopy-number rDNA: " + str(copyR)+ str("\n"))
#sys.stderr.write("[result:copy]\tCopy-number 18s: " + str(copy18)+ str("\n"))
#sys.stderr.write("[result:copy]\tCopy-number 5.8s: " + str(copy58)+ str("\n"))
#sys.stderr.write("[result:copy]\tCopy-number 28s: " + str(copy28)+ str("\n"))
#sys.stderr.write("[result:copy]\tPloidy-corrected copy-number rDNA: " + str(copyR*(float(2)/float(args['ploidy'])))+ str("\n"))
#sys.stderr.write("[result:copy]\tPloidy-corrected copy-number 18s: " + str(copy18*(float(2)/float(args['ploidy'])))+ str("\n"))
#sys.stderr.write("[result:copy]\tPloidy-corrected copy-number 5.8s: " + str(copy58*(float(2)/float(args['ploidy'])))+ str("\n"))
#sys.stderr.write("[result:copy]\tPloidy-corrected copy-number 28s: " + str(copy28*(float(2)/float(args['ploidy'])))+ str("\n"))

sys.stderr.write("[result:copy]\tAverage coverage background: " + str(meanCovHG) + str("\n"))
sys.stderr.write("[result:copy]\tAverage coverage RDNA: " + str(meanCovRDNA) + str("\n"))
sys.stderr.write("[result:copy]\tAverage coverage 18S: " + str(meanCov18) + str("\n"))
sys.stderr.write("[result:copy]\tAverage coverage 5.8S: " + str(meanCov58) + str("\n"))
sys.stderr.write("[result:copy]\tAverage coverage 28S: " + str(meanCov28) + str("\n"))

sys.stderr.write("[result:copy]\tCopy-number rDNA: " + str(copyR) + str("\n"))
sys.stderr.write("[result:copy]\tCopy-number 18s: " + str(copy18) + str("\n"))
sys.stderr.write("[result:copy]\tCopy-number 5.8s: " + str(copy58) + str("\n"))
sys.stderr.write("[result:copy]\tCopy-number 28s: " + str(copy28) + str("\n"))
#sys.stderr.write("[result:copy]\tWeighted copy-number rDNA: " + str((WcopyR*(float(2)/float(args['ploidy'])))/float(2))+ str("\n"))

resultsFileName = str(Odir + "/" + str(ID) + "_results")
out = open(resultsFileName, 'w')
outline = str("General==========\n" + "Job-ID\t" + str(ID) +"\n"+"Bam-file\t" + str(BAM)+"\n"+"Window\t" + str(args['window'])+"\n"+"Results==========\n"+"Copy-number\trDNA\t" + str(copyR) + str("\n")+"Copy-number\t18s\t" + str(copy18) + str("\n")+"Copy-number\t5.8s\t" + str(copy58)+ str("\n") + "Copy-number\t28s\t" + str(copy28)+ str("\n"))
out.write(outline)
out.close()
'''
print to stdout: avgcount(control) acgcount(rdna) avgcount(18s) avgcount(5.8s) avgcount(28s) \
copy(rdna) copy(18s) copy(5.8s) copy(28s) LOESScopy(rdna) LOESScopy(18s) LOESScopy(5.8s) LOESScopy(28s) \
diploidLOESScopy(rdna) diploidLOESScopy(18s) diploidLOESScopy(5.8s) diploidLOESScopy(28s)

print to file: LOESScoverage per window (bed)
'''
cmdline(str("rm " + bedTEMP))
cmdline(str("rm ") + windowedGCBedName)
elapsed_time = time.time() - start_time
sys.stderr.write("[info:repper]\tAll done: job " + ID + str(".\n"))
sys.stderr.write("[info:repper]\tResults are found in " + Odir + str(".\n"))
sys.stderr.write("[info:repper]\tElapsed time: "  + str(elapsed_time/60) + str(" minutes.\n"))

