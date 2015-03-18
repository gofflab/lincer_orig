#!/usr/bin/env python
'''
Created on Aug 27, 2010

@author: lgoff
'''

############
#Imports
############
import newGTF
import intervallib
import dbConn
import bisect
import sys,getopt
from misc import rstrips
import copy

############
#Constants
############
overlapThreshold = 0.20
extensionLength = 500 #grow 5'end of lincRNA by this many bases to test for Bidirectional promoter
strandLookup = {'+':'-','-':'+'}

help_message = '''
Created on Aug 27, 2010
@author: lgoff

Usage: python lincName.py [options] <gtfFile.gtf>

Options:
    -g | --genome  [Default : hg19]   Determines what build of the genome is used to fetch RefSeq transcripts
                    around which lincNames are chosen.
                    
    -h | --help       Displays this helpful help screen
    
    -v                Verbose
    
    -o | --output    [Default : <gtfFile_named.gtf>] Determines output file
'''

############
#Classes
############
class Usage(Exception):
    def __init__(self, msg):
        self.msg = msg


############
#Functions
############

def test5PrimeOverlap(lincInt,geneInt):
    """May need to validate this.  I'm not sure this works when a lincRNA completely covers a PC gene on the opposite strand"""
    assert lincInt.overlaps(geneInt)
    if lincInt.strand == "+":
        if lincInt.start <= geneInt.end and lincInt.end > geneInt.end:
            return True
        else:
            return False
    elif lincInt.strand == "-":
        if geneInt.start <= lincInt.end and geneInt.end > lincInt.end:
            return True
        else:
            return False
    else:
        raise ValueError("Could not determine")  

def bpOverlap(lincInt,geneInt):
    assert lincInt.overlaps(geneInt), "%s and %s do not overlap" % (lincInt.name,geneInt.name)
    bounds = [lincInt.start,lincInt.end,geneInt.start,geneInt.end]
    bounds.sort()
    #range = bounds[3]-bounds[0]
    overlap = bounds[2]-bounds[1]
    return overlap
        
def printLincs(handle,lincs):
    for linc in lincs:
        print >>handle, linc.getGTF(),

############
#Main
############

def main(gtfFile,genome='hg19'):
    #Parse GTF File for lincs
    lincIter = newGTF.GTFGeneIterator(gtfFile,verbose=verbose)
    
    #Retrieve and index RefSeq genes
    refSeqs = dbConn.fetchRefSeqIntervalsIndexed(genome=genome,proteinCodingOnly=True,verbose=verbose)
    
    #Results container
    res = set([])
    
    #Container for gene:linc assoc.
    geneLincs = {}
        
    #Loop through lincRNAs
    for linc in lincIter:
        flag = False
        bdFlag = False #True if linc is bidirectional
        asFlag = False #True if linc is antisense
        #Convert to Interval
        interval = linc.toInterval()
        
        #Test for weird chromosome (ie. not in refSeqs.keys() )
        if not interval.chr in refSeqs.keys():
            res.add(linc)
            continue

        #Bug tracking only
        if verbose:
            sys.stderr.write(str(interval)+"\n")
        
        #Get list of gene positions that are relevant
        senseGeneStarts = [x.start for x in refSeqs[interval.chr][interval.strand]]
        senseGeneEnds = [x.end for x in refSeqs[interval.chr][interval.strand]]
    
        #Get opposite strand to test
        testStrand = strandLookup[interval.strand]
        
        #Test overlap with genes on opposite strand
        for gene in refSeqs[interval.chr][testStrand]:
            extendedInterval = copy.copy(interval)
            extendedInterval.grow5_prime(extensionLength)
            
            if extendedInterval.overlaps(gene):
                #If 5' end of linc overlaps the 5' of a coding gene on the opposite strand, 
                #by more than 0bp but less than min(BP_THRESH * length(L), BP_THRESH * length(coding gene)) 
                #THEN name linc "linc-[HUGO_GENE_NAME]-BP"
                overlap = bpOverlap(extendedInterval,gene)
                fivePrime = test5PrimeOverlap(extendedInterval,gene)
                cutoff = min(len(extendedInterval)*overlapThreshold,gene.intervalLen()*overlapThreshold)
                if fivePrime and overlap <= cutoff:
                    linc.propogateLincName("linc-%s-BP" % gene.name)
                    linc.addAttribute("bidirectional_prom",gene.name)
                    res.add(linc)
                    flag = True
                    bdFlag = True
                    #break
                    continue
                
                #TODO FIX this so that ANY overlap that is not a BP becomes and -AS
                if not bdFlag:
                    linc.propogateLincName("linc-%s-AS" % gene.name)
                linc.addAttribute("antisense",gene.name)
                res.add(linc)
                flag = True
                asFlag = True
                break
        #ELSE find the closest coding gene on the same strand as the L, starting from the 3' end of the linc.
        #Suppose its HUGO name is NCG1.Add L to a list of lincs to be named after NCG1.
        if not flag:
            if interval.strand == "+":
                nearestGeneIdx = bisect.bisect(senseGeneStarts,interval.end) #choose most adjacent gene 3' to lincRNA
            elif interval.strand == "-":
                nearestGeneIdx = bisect.bisect(senseGeneEnds,interval.start)-1
            try:
                nearestGene = refSeqs[interval.chr][interval.strand][nearestGeneIdx]
            except IndexError:
                #If I cannot find the nearestGene (e.g. end of chromosome or something, just push linc to results
                #and deal with them later. (for now)
                
                #print nearestGeneIdx
                #print interval.toBed()
                res.add(linc)
                continue
            geneLincs.setdefault(nearestGene.name,[]).append(linc)
        
    #Evaluate container for linc:gene assocs
    """
    FOREACH coding gene G in the table above:
    IF there's only one linc to be named after G THEN
        name that linc "linc-G"
    ELSE
        sort the list of lincs by proximity to G, with the closest linc at the front of the list
        FOR i = 1 to #number of lincs named after G
            name linc i "linc-G-i"
    """
    for k,v in geneLincs.iteritems():
        if len(v) == 1:
            v[0].propogateLincName("linc-%s" % (k))
            res.add(v[0])
        elif len(v) >1:
            if v[0].strand == "+":
                v.sort(reverse=True)
            elif v[0].strand == "-":
                v.sort()
            for i in xrange(len(v)):
                v[i].propogateLincName("linc-%s-%d" % (k,i+1))
                res.add(v[i])
    return res

############
#Tests
############
def test():
    fname = '/seq/rinnscratch/cole/ftp/assemblies/linc_catalog.gtf'
    outHandle = open('/seq/rinnscratch/cole/ftp/assemblies/linc_catalog_named.gtf','w')
    verbose=True
    lincs = main(fname)
    printLincs(outHandle,lincs)
    sys.stderr.write("Done!"+"\n")
    return



############
#Orders
############
if __name__=="__main__":
    #test()
    argv = sys.argv
    #default settings
    genome = "hg19"
    verbose = False
    outFile = None
    try:
        try:
            opts,args = getopt.getopt(argv[1:],"hg:o:v",["help","genome","output"])
        except getopt.error,msg:
            raise Usage(msg)
        
        #option processing
        for option,value in opts:
            if option in ("-g","--genome"):
                genome = value
            if option in ("-h","--help"):
                raise Usage(help_message)
            if option == "-v":
                verbose = True
            if option in ("-o","--output"):
                outFile = value
        
        #debugging
        #print opts
        #print args
        
        try:        
            assert len(args)==1
            gtfFile = args[0]
        except:
            raise Usage(help_message)
        baseName = rstrips(gtfFile,".gtf")
        if verbose:
            sys.stderr.write("Naming lincs in file %s using RefSeq transcripts in genome %s.\n" % (gtfFile,genome))
        lincs = main(gtfFile,genome=genome)
        if outFile == None:
            outFile = (baseName+"_named.gtf")
        if verbose:
            sys.stderr.write("Writing output to %s.\n" % outFile)
        outHandle = open(outFile,'w')
        printLincs(outHandle,lincs)
        if verbose:
            sys.stderr.write("Done!\n")
    except Usage, err:
        print >>sys.stderr, sys.argv[0].split("/")[-1] + ": " + str(err.msg)
        sys.exit()
    