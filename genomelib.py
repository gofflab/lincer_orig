'''
Created on Aug 28, 2010

This is a port of the genome.py module from seqtools (it is a work in progress)

@author: lgoff
'''
############
#Imports
############
import sequencelib
import random
#######
#Constants
#######

purines=['A','G']
pyrimidines=['C','T','U']

chr_names = ['chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10',
             'chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19',
             'chr20','chr21','chr22','chrX','chrY']

genome_length = 3080419480

chr_lengths = {'chr1':247249719,
               'chr2':242951149,
               'chr3':199501827,
               'chr4':191273063,
               'chr5':180857866,
               'chr6':170899992,
               'chr7':158821424,
               'chr8':146274826,
               'chr9':140273252,
               'chr10':135374737,
               'chr11':134452384,
               'chr12':132349534,
               'chr13':114142980,
               'chr14':106368585,
               'chr15':100338915,
               'chr16':88827254,
               'chr17':78774742,
               'chr18':76117153,
               'chr19':63811651,
               'chr20':62435964,
               'chr21':46944323,
               'chr22':49691432,
               'chrX':154913754,
               'chrY':57772954
               }

genbases = {'A': 843953565, 'C': 584268578, 'T': 845168978, 'G': 584621685, 'N': 222406671}
genfreqs = {'A': 0.27397358394837834, 'C': 0.18967175795161509, 'T': 0.27436814482162669, 'G': 0.18978638746954035, 'N': 0.072200124834946186}
genome_build = 'hg18'
genome_dir = '/seq/compbio-hp/lgoff/genomes/'+genome_build
genome_file = genome_build+".fa"
hg19_genome_file = '/fg/compbio-t/lgoff/magda/references/human/genome/hg19/hg19.fa'
hg18_genome_file = '/fg/compbio-t/lgoff/magda/references/human/genome/hg18/hg18.fa'

rmgenome_dir = "/seq/compbio-hp/lgoff/smallRNAs/genomes/human_repeatmasked/"

bed_fields = ['chr','start','end','label','score','strand']

mammals_alignments_dir = '/ahg/scr3/mammals/ucsc/multiz44way/'

#######
#Functions
#######
def fetch_genbases(genhandle,genbases={}):
    bases = ['A','T','G','C','N']
    geniter = sequencelib.FastaIterator(genhandle)
    for genseq in geniter:
        print genseq['name']
        seq = genseq['sequence'].upper()
        for b in bases:
            genbases[b] = seq.count(b) + genbases.get(b,0)
    return genbases

def fetch_genome_freqs():
    """Specifically returns a dictionary containing frequencies of every 7mer in hg18"""
    freqfile = '/seq/compbio-hp/lgoff/smallRNAs/genomes/human/hg18/hg18_7mer_frequencies.txt'
    freqhandle = open(freqfile,'r')
    freqs = {}
    for line in freqhandle:
        vals = line.rstrip().split()
        freqs[vals[0]] = float(vals[1])
    return freqs


def random_region(n,m=1):
    '''Generate a random region of max length "n" and min length "m" (default m=1).'''
    c = random.choice(chr_names)
    strand= random.choice(["+","-"])
    start = random.randint(1,chr_lengths[c])
    end = start+random.randint(m,n)
    return c, start, end, strand

def isMasked(s):
    maskedChars='actgnN'
    for c in s:
        if c in maskedChars:
            return True
    return False


#######################
#pygr specific
#######################
#SeqPath = pygr.Data.Bio.Seq.Genome.HUMAN.hg18

def pygrConnect(genome="hg18",useWorldbase = True):
    if genome in ["hg18","hg19"]: useWorldbase = False
    if useWorldbase:
        from pygr import worldbase
        if genome == "hg18":
            res=worldbase.Bio.Seq.Genome.HUMAN.hg18()
        elif genome == "hg19":
            res=worldbase.Bio.Seq.Genome.HUMAN.hg19()
        elif genome == "mm9":
            res=worldbase.Bio.Seq.Genome.MOUSE.mm9()
        elif genome == "mm8":
            res=worldbase.Bio.Seq.Genome.MOUSE.mm8()
        else:
            raise AssertionError ("No genome by that name in worldbase. (I think)...")
    else:
        from pygr import seqdb
        if genome == "hg18":
            res = seqdb.SequenceFileDB('/fg/compbio-t/lgoff/magda/references/human/genome/hg18/hg18.fa')
        elif genome == "hg19":
            res = seqdb.SequenceFileDB('/n/rinn_data1/indexes/human/hg19/hg19.fa')
        elif genome == "mm9":
            res = seqdb.SequenceFileDB("/n/rinn_data1/indexes/igenomes/Mus_musculus/UCSC/mm9/Sequence/Chromosomes/mm9.fa")
        else:
            raise AssertionError ("I'm not sure how to handle that genome build yet...sorry.")
    return res

def pygrConnectBroad(genome="hg18",useWorldbase=True):
    if genome in ["hg18","hg19"]: useWorldbase = False
    if useWorldbase:
        from pygr import worldbase
        if genome == "hg18":
            res=worldbase.Bio.Seq.Genome.HUMAN.hg18()
        elif genome == "hg19":
            res=worldbase.Bio.Seq.Genome.HUMAN.hg19()
        elif genome == "mm9":
            res=worldbase.Bio.Seq.Genome.MOUSE.mm9()
        elif genome == "mm8":
            res=worldbase.Bio.Seq.Genome.MOUSE.mm8()
        else:
            raise AssertionError ("No genome by that name in worldbase. (I think)...")
    else:
        from pygr import seqdb
        if genome == "hg18":
            res = seqdb.SequenceFileDB('/fg/compbio-t/lgoff/magda/references/human/genome/hg18/hg18.fa')
        elif genome == "hg19":
            res = seqdb.SequenceFileDB('/fg/compbio-t/lgoff/magda/references/human/genome/hg19/hg19.fa')
        else:
            raise AssertionError ("I'm not sure how to handle that genome build yet...sorry.")
    return res

def fetchSequence(chrom,start,end,strand):
    hg18=pygrConnect()
    start,end=int(start),int(end)
    seq=hg18[chrom][start:end]
    if strand == "-":
        seq=-seq
    return seq
