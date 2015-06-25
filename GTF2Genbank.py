#!/usr/bin/env python
# encoding: utf-8
"""
GTF2Genbank.py

Created by Loyal Goff on 2010-11-08.
Copyright (c) 2010 __MyCompanyName__. All rights reserved.
"""

import sys
import getopt
import newGTF,genomelib,dbConn
import datetime,string
from misc import rstrips
import time


help_message = '''
The help message goes here.
'''


class Usage(Exception):
    def __init__(self, msg):
        self.msg = msg

def printHeader(outHandle,gtf_gene_container,date,verbose=False):
    try:
        lincName = gtf_gene_container.transcripts[0].features[0].attributes['linc_name']
    except:
        lincName = gtf_gene_container.transcripts[0].features[0].attributes['gene_id']

    print >>outHandle, "LOCUS\t%s\t%d bp\tDNA\tPRI\t%s" % (lincName,len(gtf_gene_container.sequence),date)
    print >>outHandle, "DEFINITION\tlincRNA genomic sequence for %s from cufflinks assembly." % (lincName)
    print >>outHandle, "ACCESSION\t%s" % (gtf_gene_container.geneId)
    if verbose:
        sys.stderr.write("\tHeader\t")
    return

def printSequence(outHandle,gtf_gene_container,verbose=False):
    blockSize = 10    
    seq = gtf_gene_container.sequence
    #sys.stderr.write("%d" % len(seq))
    blocks = []
    printStringElements = []
    for k in xrange(0,len(seq),blockSize):
        blocks.append(seq[k:k+blockSize])
    #sys.stderr.write("Got Blocks\t")
    for i in xrange(len(blocks)):
        if i%6 == 0:
            #sys.stderr.write("%d" % i)
            printStringElements.append("\n%s " % string.rjust(str((i*10+1)),10))
        printStringElements.append("%s " % blocks[i])
    printString = ''.join(printStringElements)
    print >>outHandle, "ORIGIN"
    print >>outHandle, printString
    if verbose:
        sys.stderr.write("Sequence\n")
    return

def printFeatures(outHandle,gtf_gene_container,genome='hg19',verbose=False):
    try:
        lincName = gtf_gene_container.transcripts[0].features[0].attributes['linc_name']
    except:
        lincName = gtf_gene_container.transcripts[0].features[0].attributes['gene_id']
    print >>outHandle, "FEATURES\t\tLocation/Qualifiers"
    print >>outHandle, '\tsource\t1..%d\n\t\t/organism="%s"\n\t\t/build="%s"\n\t\t/db_xref="taxon:9606"\n\t\t/chromosome="%s"' % (len(gtf_gene_container.sequence),genome,genome,gtf_gene_container.contig)
    if gtf_gene_container.strand=="+":
        print >>outHandle, '\tgene\t1..%d\n\t\t/gene="%s"\n\t\t/xloc="%s"' % (len(gtf_gene_container.sequence),lincName,gtf_gene_container.geneId)
    elif gtf_gene_container.strand=="-":
        print >>outHandle, '\tgene\tcomplement(1..%d)\n\t\t/gene="%s"\n\t\t/xloc="%s"' % (len(gtf_gene_container.sequence),lincName,gtf_gene_container.geneId)
    for transcript in gtf_gene_container.transcripts:
        interval = transcript.toSplicedInterval()
        joinString = 'join(%s)' % ",".join([str(interval.exonStarts[i]-gtf_gene_container.start+2)+".."+str(interval.exonEnds[i]-gtf_gene_container.start+1) for i in xrange(len(interval.exonStarts))])
        if interval.strand == "-":
            joinString = 'complement('+joinString+')'
        print >>outHandle, '\tncRNA\t%s' % joinString
        print >>outHandle, '\t\t/gene="%s"' % interval.name
        """
        try: print >>outHandle, '\t\t/class_code="%s"' % transcript.features[0].attributes['class_code']
        except: pass
        try: print >>outHandle, '\t\t/linc_name="%s"' % transcript.features[0].attributes['linc_name']
        except: pass
        try: print >>outHandle, '\t\t/nearest_ref="%s"' % transcript.features[0].attributes['nearest_ref']
        except: pass
        try: print >>outHandle, '\t\t/oId="%s"' % transcript.features[0].attributes['nearest_ref']
        except: pass
        try: print >>outHandle, '\t\t/gene_id="%s"' % transcript.features[0].attributes['gene_id']
        except: pass
        """
        #OR Try This
        for k in transcript.features[0].attributes.keys():
            if not k.startswith('exon'):
                print >>outHandle, '\t\t/%s="%s"' % (k,transcript.features[0].attributes[k])
        
    if verbose:
        sys.stderr.write("Features\t")
    return

def printRmskFeatures(outHandle,interval,cursor=None):
    rmsk = dbConn.findRepeatOverlap(interval,cursor=cursor)
    if rmsk:
        for r in rmsk:
            start= int(r['genoStart'])-interval.start+1
            end = int(r['genoEnd'])-interval.start
            strand = r['strand']
            coordString = ''
            if start < 1:
                coordString += '<1..'
            else:
                coordString += '%d..' % (start)
            if end > len(interval)+1:
                coordString += '%d>' % (len(interval))
            else:
                coordString += '%d' % (end)
            if strand == "-":
                coordString = 'complement(%s)' % (coordString)
            print >>outHandle, "\trepeat_region\t%s" % (coordString)
            print >>outHandle, '\t\t/rpt_type="%s"\n\t\t/rpt_family="%s"\n\t\t/standard_name="%s"' % (r['repClass'],r['repFamily'],r['repName'])
    return

def printUCSCFeatures(outHandle,interval,cursor=None):
    ucsc = dbConn.findUCSCOverlap(interval,cursor=cursor)
    if ucsc:
        #print len(ucsc)
        for u in ucsc:
            start = int(u['txStart'])-interval.start+1
            end = int(u['txEnd'])-interval.start
            strand = u['strand']
            exonStarts = [int(x)-interval.start+1 for x in u['exonStarts'].split(",")[:-1]]
            exonEnds = [int(x)-interval.start for x in u['exonEnds'].split(",")[:-1]]
            #DEBUGGING
#            print exonStarts
#            print u['exonStarts']
#            print exonEnds
#            print u['exonEnds']
            
            joinString = ''
            for i in xrange(len(exonStarts)):
                tmpJoinString = ''
                #TODO: Check to make sure cases cannot overlap
                if exonStarts[i] < 1 and exonEnds[i] < 1:
                    continue
                
                elif exonStarts[i] > len(interval)+1 and exonEnds[i] >len(interval)+1:
                    continue
                
                elif exonStarts[i] < 1 and exonEnds[i] <= len(interval)+1:
                    tmpJoinString += '<1..%d' % (exonEnds[i])
                    
                elif exonStarts[i] >=1 and exonEnds[i] <= len(interval)+1:
                    tmpJoinString += '%d..%d' % (exonStarts[i],exonEnds[i])
                    
                elif exonStarts[i] >=1 and exonEnds[i] > len(interval)+1:
                    tmpJoinString += '%d..%d>' % (exonStarts[i],len(interval))
                
                if strand == "-":
                    tmpJoinString = 'complement(%s)' % (tmpJoinString)
                joinString += tmpJoinString+","
            if not joinString == '':
                joinString = 'join(%s)' % (joinString.rstrip(","))
                print >>outHandle, '\tmRNA\t%s' % (joinString)
                print >>outHandle, '\t\t/gene="%s"' % (u['name']) 
    return

#####
#Utilities
#####
def fetch_bases(seq,basecount={}):
    bases = ['A','T','G','C','N']
    seq = seq['sequence'].upper()
    for b in bases:
        basecount[b] = seq.count(b) + basecount.get(b,0)
    return basecount

########
#main
########

def main(argv=None):
    if argv is None:
        argv = sys.argv
    try:
        try:
            opts, args = getopt.getopt(argv[1:], "hg:o:s:v", ["help", "genome=","output=","system="])
        except getopt.error, msg:
            raise Usage(msg)
        #Defaults
        date = datetime.date.today()
        genome = 'hg19'
        outFile = None
        verbose = False
        system = ''
        # option processing
        for option, value in opts:
            if option == "-v":
                verbose = True
            if option in ("-h", "--help"):
                raise Usage(help_message)
            if option in ("-o", "--output"):
                outFile = value
            if option in ("-g", "--genome"):
                genome = value
            if option in ("-s", "--system"):
                system = value
        
        assert len(args)==1
        fname = args[0]
        if outFile == None:
            outFile = rstrips(fname,".gtf")+".gbk"
        outHandle = open(outFile,'w')
        iter = newGTF.GTFGeneIterator2(fname)
        genes = []
        if verbose:
            sys.stderr.write("Fetching Genes from gtf file.\n")
        for i in iter:
            genes.append(i)
        geneCount = 0
        if system == "Broad":
            pygrConnection = genomelib.pygrConnectBroad(genome,useWorldbase=False)
            gbdbCursor = dbConn.gbdbConnect(gbdbname=genome)
        else:
            pygrConnection = genomelib.pygrConnect(genome,useWorldbase=False)
            gbdbCursor = dbConn.gbdbConnect(gbdbname=genome)
        for gene in genes:
            interval = gene.toInterval()
            interval.genome = genome
            geneCount +=1
            sys.stderr.write("%d " % geneCount)
            if verbose:
                sys.stderr.write("%s\n" % gene.geneId)
            #gene.fetchSequence(genome)
            gene.fetchSequence(genome,pygrConnection)
            if verbose:
                sys.stderr.write("Got Sequence\t")
            seqLen = len(gene.sequence)
            printHeader(outHandle,gene,date,verbose=verbose)
            printFeatures(outHandle,gene,genome=genome,verbose=verbose)
            printRmskFeatures(outHandle,interval,cursor=gbdbCursor)
            printUCSCFeatures(outHandle,interval,cursor=gbdbCursor)
            printSequence(outHandle,gene,verbose=verbose)
            print >>outHandle, "\n//\n"    
        
        outHandle.close()
    except Usage, err:
        print >> sys.stderr, sys.argv[0].split("/")[-1] + ": " + str(err.msg)
        print >> sys.stderr, "\t for help use --help"
        return 2

if __name__ == "__main__":
    sys.exit(main())
