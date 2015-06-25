#!/usr/bin/env python
import MySQLdb,sys,time
import intervallib
import genomelib
import sequencelib

###################
#
#Connect to Broad MySQL Database
#
###################
def broadConnect():
    host="mysql.broadinstitute.org"
    user="lgoff"
    password="nextgen"
    db="lgoff_nextgen"
    broadDb=MySQLdb.connect(host=host,user=user,db=db,passwd=password)
    return broadDb.cursor(MySQLdb.cursors.DictCursor)
    
###################
#
#Connection to UCSC Genome Browser MySQL Database
#
###################
def gbdbConnect(gbdbname = "hg19"):
    gbHost = "genome-mysql.cse.ucsc.edu"
    gbUser = "genome"
    gbdb = MySQLdb.connect(host=gbHost,user=gbUser,db=gbdbname)
    return gbdb.cursor(MySQLdb.cursors.DictCursor)

###################
#
#Connection to Valor local UCSC Genome Browser MySQL Database
#
###################
def valorGbdbConnect(gbdbname='hg19'):
    gbHost = 'localhost'
    gbUser = 'root'
    gbPass = 'rinnlab1'
    gbdb = MySQLdb.connect(host=gbHost,user=gbUser,passwd=gbPass,db=gbdbname)
    return gbdb.cursor(MySQLdb.cursors.DictCursor)

###################
#
#Connection to Ensembl MySQL Database
#
####################
def ensemblConnect():
    ensemblHost = "ensembldb.ensembl.org"
    ensemblUser = "anonymous"
    ensembldbname = "homo_sapiens_core_47_36i"
    ensembldb = MySQLdb.connect(host=ensemblHost,user=ensemblUser,db=ensembldbname)
    return ensembldb.cursor(MySQLdb.cursors.DictCursor)

####################
#
#Operations on UCSC genome browser data
#
####################
def fetchRefSeq(genome = 'hg18',lookupval = 'name'):
    """Returns a dictionary of RefSeq genes (by chromosome and strand with 'name' parameter as key) from UCSC genome browser (equivalent to RefSeq ID)"""
    cursor=gbdbConnect(gbdbname=genome)
    select="SELECT * FROM refGene"
    cursor.execute(select)
    rows=cursor.fetchall()
    output={}
    for chr in genomelib.chr_names:
        output[chr]={}
        output[chr]['+']={}
        output[chr]['-']={}
    for row in rows:
        if row['chrom'] in genomelib.chr_names:
            output[row['chrom']][row['strand']][row[lookupval]]=row
    return output 

def fetchRefSeqIntervals(genome = 'hg18'):
    cursor = gbdbConnect(gbdbname=genome)
    select = "SELECT * from refGene"
    cursor.execute(select)
    rows = cursor.fetchall()
    output = {}
    for row in rows:
        exonStarts = map(int,row['exonStarts'].rstrip().split(","))
        exonEnds = map(int,row['exonEnds'].rstrip().split(","))
        start = int(row['txStart'])
        exonOffsets = [x-start for x in exonStarts]
        exonLengths = []
        for i in len(exonStarts):
            exonLengths.append(exonEnds-exonStarts+1)
        output[row['name']] = intervallib.SplicedInterval(row['chrom'],row['txStart'],row['txEnd'],row['strand'],",".join([str(x) for x in exonLengths]),",".join([str(x) for x in exonOffsets]),name=row['name2'])
    return output

def fetchRefSeqIntervalsIndexed(genome='hg18',proteinCodingOnly=False,verbose=False):
    """
    Returns a dictionary of RefSeq SplicedIntervals (by chromosome and strand) from UCSC table browser.
    Indexed lists are sorted prior to return for easy search
    Same as fetchRefSeqIntervals but indexed by chrom and strand
    """
    cursor=gbdbConnect(gbdbname=genome)
    select="SELECT * FROM refGene"
    if verbose:
        sys.stderr.write("Fetching RefSeq Sequences...\n")
    cursor.execute(select)
    rows=cursor.fetchall()
    output={}
    for chr in genomelib.chr_names:
        output[chr]={}
        output[chr]['+']=[]
        output[chr]['-']=[]
    if verbose:
        sys.stderr.write("Creating index by chr and strand...\n")
    
    for row in rows:
        if proteinCodingOnly and not row['name'].startswith('NM'):
            continue
        try:
            exonStarts = map(int,row['exonStarts'].rstrip().split(",")[:-1])
            exonEnds = map(int,row['exonEnds'].rstrip().split(",")[:-1])
        except:
            print "\t".join(["%s:%s" % (k,v) for k,v in row.iteritems()])
        start = int(row['txStart'])
        exonOffsets = [x-start for x in exonStarts]
        exonLengths = []
        for i in xrange(len(exonStarts)):
            exonLengths.append(exonEnds[i]-exonStarts[i]+1)
        if row['chrom'] in genomelib.chr_names:
            output[row['chrom']][row['strand']].append(intervallib.SplicedInterval(row['chrom'],row['txStart'],row['txEnd'],row['strand'],",".join([str(x) for x in exonLengths]),",".join([str(x) for x in exonOffsets]),name=row['name2']))
    
    #Sort 
    if verbose:
        sys.stderr.write("Sorting:\n")
    tstart = time.time()
    for key in output.keys():
        if verbose:
            sys.stderr.write("\t%s\t" % key)
        output[key]['+'].sort()
        output[key]['-'].sort()
        tend = time.time()
        if verbose:
            sys.stderr.write('%0.2f sec\n' % (tend-tstart))
        tstart = time.time()
    return output

def getIntervalFromRefSeq(lookupval,genome='hg18',lookupkey= 'name2',verbose=False):
    cursor = gbdbConnect(gbdbname=genome)
    select = """SELECT * FROM refGene WHERE %s = '%s'""" % (lookupkey,lookupval)
    if verbose:
        sys.stderr.write("Query: "+select+"\nFetching RefSeq Record(s)\n")
    cursor.execute(select)
    rows=cursor.fetchall()
    if verbose:
        sys.stderr.write("%d Rows returned...\n" % len(rows))
    output = []
    for row in rows:
        try: 
            exonStarts = map(int,row['exonStarts'].rstrip().split(",")[:-1])
            exonEnds = map(int,row['exonEnds'].rstrip().split(",")[:-1])
        except:
            print "\t".join(["%s:%s" % (k,v) for k,v in row.iteritems()])
        start = int(row['txStart'])
        exonOffsets = [x-start for x in exonStarts]
        exonLengths = []
        for i in xrange(len(exonStarts)):
            exonLengths.append(exonEnds[i]-exonStarts[i]+1)
        output.append(intervallib.SplicedInterval(row['chrom'],row['txStart'],row['txEnd'],row['strand'],",".join([str(x) for x in exonLengths]),",".join([str(x) for x in exonOffsets]),name=row['name2']))
    return output

def getIntervalFromAll_mRNA(lookupval,genome='hg18',lookupkey='qName',verbose=False):
    cursor = gbdbConnect(gbdbname=genome)
    select = """SELECT * FROM all_mrna WHERE %s = '%s'""" % (lookupkey,lookupval)
    if verbose:
        sys.stderr.write("Query: "+select+"\nFetching all_mrna Record(s)\n")
    cursor.execute(select)
    rows=cursor.fetchall()
    if verbose:
        sys.stderr.write("%d Rows returned...\n" % len(rows))
    output = []
    for row in rows:
        try:
            exonStarts = map(int,row['tStarts'].rstrip().split(",")[:-1])
            blockSizes = map(int,row['blockSizes'].rstrip().split(",")[:-1])
            exonEnds = [exonStarts[i]+blockSizes[i] for i in xrange(len(exonStarts))]
        except:
            print "\t".join(["%s:%s" % (k,v) for k,v in row.iteritems()])
        start = int(row['tStart'])
        exonOffsets = [x-start for x in exonStarts]
        exonLengths = [exonEnds[i]-exonStarts[i]+1 for i in xrange(len(exonStarts))]
        output.append(intervallib.SplicedInterval(row['tName'],start,int(row['tEnd']),row['strand'],",".join([str(x) for x in exonLengths]),",".join([str(x) for x in exonOffsets]),name=row['qName']))
    return output

def refseqTSS():
    """Uses fetchRefSeq to retrieve current RefSeq Sequences and then returns a sorted list of tuples (as value of chr.strand dictionaries) containing ('refSeqID','chr','tss','orientation')"""
    refSeqs=fetchRefSeq()
    output={}
    for chr in genomelib.chr_names:
        output[chr]=[]
        for strand in ['+','-']:
            for k in refSeqs[chr][strand]:
                v=refSeqs[chr][strand][k]
                if v['strand'] == "+":
                    tss=v['txStart']
                elif v['strand'] == "-":
                    tss=v['txEnd']
                tssInfo=(v['name'],v['chrom'],int(tss),v['strand'])
                output[chr].append(tssInfo)
            output[chr].sort(lambda x,y:cmp(x[2],y[2]))
    return output

def fetchwgRNA():
    cursor=gbdbConnect()
    select="SELECT * FROM wgRna"
    cursor.execute(select)
    rows=cursor.fetchall()
    output={}
    for chr in genomelib.chr_names:
        output[chr]={}
        output[chr]['+']={}
        output[chr]['-']={}
    for row in rows:
        if row['chrom'] in genomelib.chr_names:
            output[row['chrom']][row['strand']][row['name']]=row
    return output
#Tests for known annotation
def hostRefSeq(chr,start,end,strand):
    """
    Checks to see if interval is within a host RefSeq gene (does not test strand!!).  If no, returns False.  
    If yes, returns a list of dictionaries for each host RefSeq gene.  Keys are consistent with field names 
    from UCSC table refGene.
    """
    cursor=gbdbConnect()
    selSQL="SELECT * from refGene WHERE chrom='%s' AND txStart<='%d' AND txEnd>='%d'" % (chr,int(start),int(end))
    cursor.execute(selSQL)
    rows=cursor.fetchall()
    results=[]
    if cursor.rowcount==0:
        return False
    else:
        for row in rows:
            results.append(row)
        return results

def testCpG(chr,start,end):
    cursor=gbdbConnect()
    selSQL="SELECT * from cpgIslandExt WHERE chrom='%s' AND chromStart<='%d' AND chromEnd>='%d'" % (chr,int(start),int(end))
    cursor.execute(selSQL)
    if cursor.rowcount==0:
        return False
    else:
        return cursor.fetchone()

def testwgRNA(chr,start,end,strand):
    """
    Checks to see if interval is entirely within a known wgRNA gene (including miRNA). Does consider strand!!!
    If no flanking host wgRNA, returns False. If yes, returns a list of dictionaries for each host wgRNA gene.
    Keys are consistent with field names from UCSC table wgRNA.
    """
    cursor=gbdbConnect()
    selSQL="SELECT * from wgRna WHERE chrom='%s' AND strand='%s' AND chromStart<='%d' AND chromEnd>='%d'" % (chr,strand,int(start),int(end))
    cursor.execute(selSQL)
    rows=cursor.fetchall()
    results=[]
    if cursor.rowcount==0:
        return False
    else:
        for row in rows:
            results.append(row)
        return results

def hostmRNA(chr,start,end,strand):
    cursor=gbdbConnect()
    selSQL="SELECT * from %s_mrna WHERE tName='%s' AND tStart<='%d' AND tEnd>='%d'" % (chr,chr,int(start),int(end))
    cursor.execute(selSQL)
    rows=cursor.fetchall()
    results=[]
    if cursor.rowcount==0:
        return False
    else:
        for row in rows:
            results.append(row)
        return results

def fetchLincRNA(fname="/seq/compbio/lgoff/lincRNAs/hg18_lincRNA_Guttman.bed"):
    handle=open(fname,'r')
    lincs={}
    for chr in genomelib.chr_names:
        lincs[chr]=[]
    for line in handle:
        if line.startswith("#"):continue
        fields=['chr','start','end']
        vals=line.rstrip().split("\t")
        d=dict(zip(fields,vals))
        d['start'],d['end']=int(d['start']),int(d['end'])
        lincs[d['chr']].append(d)
    return lincs

def fetchmiRNASeeds(fname="/seq/compbio/lgoff/smallRNAs/genomes/human/microRNA/mature.fa",species = 'hsa'):
    handle = open(fname,'r')
    seeds = {}
    iter = sequencelib.FastaIterator(handle)
    for i in iter:
        if i.name.startswith(species):
            seeds[i.sequence[1:8]] = i.name.split()[0]
    return seeds

#############
#Added for lincRNA pipeline (only works on valor)
############

def findRepeatOverlap(interval,cursor=None):
    if cursor == None:
        cursor = valorGbdbConnect(interval.genome)
    if interval.genome == "mm9":
        selSQL = "SELECT * from %s_rmsk WHERE genoName = '%s' AND (genoStart >= '%d' OR genoEnd >= '%d') AND (genoStart <= '%d' OR genoEnd <= '%d')" % (interval.chr, interval.chr,interval.start,interval.start,interval.end,interval.end)
    else:
        selSQL = "SELECT * from rmsk WHERE genoName = '%s' AND (genoStart >= '%d' OR genoEnd >= '%d') AND (genoStart <= '%d' OR genoEnd <= '%d')" % (interval.chr,interval.start,interval.start,interval.end,interval.end)
    cursor.execute(selSQL)
    rows = cursor.fetchall()
    results=[]
    if cursor.rowcount==0:
        return False
    else:
        for row in rows:
            results.append(row)
        return results
    
def findUCSCOverlap(interval,cursor=None):
    if cursor == None:
        cursor = valorGbdbConnect(interval.genome)
    selSQL = "SELECT * from knownGene kg LEFT JOIN knownToRefSeq krs ON kg.name = krs.name WHERE kg.chrom = '%s' AND (kg.txStart >= '%d' OR kg.txEnd >= '%d') AND (kg.txStart <= '%d' OR kg.txEnd <= '%d')" % (interval.chr,interval.start,interval.start,interval.end,interval.end)
    cursor.execute(selSQL)
    rows = cursor.fetchall()
    results = []
    if cursor.rowcount == 0:
        return False
    else:
        for row in rows:
            results.append(row)
        return results
