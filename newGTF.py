'''
Created on Aug 31, 2010

All of this is very fragile and is
absolutely dependent on a unique geneId and unique transcriptId for any records...

@author: lgoff
'''
###########
#Imports
###########
import intervallib
import sys
from misc import uniqify,pp
import genomelib

#######################
#Error Handling
#######################
class Error(Exception):
    """Base class for exceptions in this module."""
    def __str__(self):
        return str(self.message)
    def _get_message(self, message): return self._message
    def _set_message(self, message): self._message = message
    message = property(_get_message, _set_message)

class ParsingError(Error):
    """
    Exception raised for errors in the input.

    Attributes:
        message -- explanation of the error
    """
    def __init__(self, message):
        self.message = message

#########################
#GTF Entry Class
#########################

class GTF_Entry:
    '''
    Holds a row's worth of GTF information.
    '''

    def __init__(self):
        '''
        Constructor
        '''
        self.contig = "."
        self.source = "."
        self.feature = "."
        self.frame = "."
        self.start = 0
        self.end = 0
        self.score = "."
        self.strand = "."
        self.attributes = {}

    def __cmp__(self,b):
        mid1 = (self.start+self.end)/2
        mid2 = (b.start+b.end)/2
        return cmp(mid1,mid2)

    def __repr__(self):
        return self.attributes['transcript_id']+":"+self.feature

    def addGTF_Entry(self,gtf_entry):
        self.contig = gtf_entry.contig
        self.source = gtf_entry.source
        self.feature = gtf_entry.feature
        self.frame = gtf_entry.frame
        self.start = int(gtf_entry.start)
        self.end = int(gtf_entry.end)
        self.score = gtf_entry.score
        self.strand = gtf_entry.strand
        self.attributes = gtf_entry.attributes

    def read(self,line):
        """
        read gff entry from line.
        <seqname> <source> <feature> <start> <end> <score> <strand> <frame> [attributes] [comments]
        """
        data = line.rstrip().split("\t")

        try:
            (self.contig, self.source, self.feature,
             self.start, self.end, self.score, self.strand,
             self.frame ) = [x.strip() for x in data[:8]]
        except ValueError:
            raise ValueError( "parsing error in line `%s`" % line )

        ## note: frame might be "."
        (self.start, self.end) = map(int, (self.start, self.end))
        try:
            self.score = float(self.score)
        except:
            pass
        #TODO: This may only be necessary when I convert to an Interval object
        #self.start -= 1

        self.parseInfo( data[8], line )

    def parseInfo(self,myAttributes,line ):
        """
        Parse attributes.
        """
        # remove comments
        myAttributes = myAttributes.split( "#" )[0]
        # separate into fields
        fields = map( lambda x: x.strip(), myAttributes.split(";")[:-1])
        self.attributes = {}

        for f in fields:
            d = map( lambda x: x.strip(), f.split(" "))
            n,v = d[0], d[1]
            if len(d) > 2: v = d[1:]
            if v[0] == '"' and v[-1] == '"':
                v = v[1:-1]
            else:
                ## try to convert to a value
                try:
                    v = float( v )
                    v = int( v )
                except ValueError:
                    pass
                except TypeError:
                    pass
            self.attributes[n] = v

    def toGTF(self):
        tmp = '%s\t%s\t%s\t%d\t%d\t%s\t%s\t%s\t' % (self.contig,self.source,self.feature,self.start,self.end,str(self.score),self.strand,self.frame)
        #Print 'gene_id' and 'transcript_id' as first and second attributes (required by GTF spec.)
        for attr in ['gene_id','transcript_id']:
            try:
                tmp += '%s "%s"; ' % (attr,self.attributes[attr])
            except:
                pass
        #Print remainder of attributes in any order.
        for k,v in self.attributes.iteritems():
            if k not in ['gene_id','transcript_id']:
                tmp += '%s "%s"; ' % (k,str(v))
        tmp += "\n"
        return tmp

############
#GTFTranscriptContainer
############
class GTFTranscriptContainer(object):
    def __init__(self):
        '''
        Constructor
        '''
        self.features = []
        self.start = -1
        self.end = -1
        self.contig = None
        self.strand = None
        self.transcriptId = ''
        self.geneId = ''

    def __len__(self):
        return self.end-self.start+1

    def __cmp__(self,b):
        mid1 = (self.start+self.end)/2
        mid2 = (b.start+b.end)/2
        return cmp(mid1,mid2)

    def addFeature(self,gtf_entry):
        if self.transcriptId == '':
            self.contig = gtf_entry.contig
            self.strand = gtf_entry.strand
            self.transcriptId = gtf_entry.attributes['transcript_id']
        assert self.transcriptId == gtf_entry.attributes['transcript_id']
        if 'gene_id' in gtf_entry.attributes.keys():
            self.geneId = gtf_entry.attributes['gene_id']
        self.features.append(gtf_entry)
        self.update()

    def update(self):
        self.start = min([x.start for x in self.features])
        self.end = max([x.end for x in self.features])

    def toSplicedInterval(self):
        transcripts = uniqify([x.attributes['transcript_id'] for x in self.features])
        if len(transcripts) > 1:
            raise ValueError ("Something is wrong, there are too many different transcript_ids")
        for t in transcripts:
            exons = [x for x in self.features if (x.feature =='exon' and x.attributes['transcript_id']==t)]
            exons.sort(key=lambda x: x.attributes['exon_number'])
            transStart = min([x.start-1 for x in exons])
            myInt = intervallib.SplicedInterval(self.contig,transStart,max([x.end for x in exons]),self.strand,",".join([str(x.end-x.start+1) for x in exons]),",".join([str(x.start-1-transStart) for x in exons]),name=t)
            return myInt


############
#Gene Container
############

class GTFGeneContainer(object):
    '''
    Container for all GTF_Entry instances with a common geneId
    Assumptions:
        - gene_id field is unique to a gene locus (ie. not shared amongst gene duplicates
        - There is no guarantee that the order of rows is preserved during reading in and returning GTF

    '''

    def __init__(self):
        '''
        Constructor
        '''
        self.features = []
        self.transcripts = []
        self.start = -1
        self.end = -1
        self.contig = None
        self.strand = None
        self.geneId = ''
        self.sequence = ''

    def __len__(self):
        return self.end-self.start+1

    def __cmp__(self,b):
        mid1 = (self.start+self.end)/2
        mid2 = (b.start+b.end)/2
        return cmp(mid1,mid2)

    def addFeature(self,gtf_entry):
        if self.geneId == '':
            self.contig = gtf_entry.contig
            self.strand = gtf_entry.strand
            self.geneId = gtf_entry.attributes['gene_id']
        assert self.geneId == gtf_entry.attributes['gene_id']
        self.features.append(gtf_entry)
        self.update()

    def addGTFTranscript(self,gtf_transcript):
        if self.geneId == '':
            self.contig = gtf_transcript.contig
            self.strand = gtf_transcript.strand
            self.geneId = gtf_transcript.geneId
        assert self.geneId == gtf_transcript.geneId and self.contig == gtf_transcript.contig and self.strand == gtf_transcript.strand
        self.transcripts.append(gtf_transcript)
        self.transcriptUpdate()

    def update(self):
        self.start = min([x.start for x in self.features])
        self.end = max([x.end for x in self.features])

    def transcriptUpdate(self):
        self.start = min([x.start for x in self.transcripts])
        self.end = max([x.end for x in self.transcripts])


    def propogateLincName(self,lincName):
        for feat in self.features:
            feat.attributes['linc_name'] = lincName
            #if not 'gene_name' in feat.attributes:
            feat.attributes['gene_name'] = lincName

    def addAttribute(self,key,value):
        for feat in self.features:
            feat.attributes[key] = value

    def geneToBed(self):
        """This does not work yet"""
        raise Error ("This method does not work yet")
        return "%s\t%d\t%d\t%s\t0\t%s\t%s\t%s" % (self.contig,self.start,self.end,self.attributes['transcript_id'],self.strand,",".join(self.exonLengths),",".join(self.exonOffsets))

    def transcriptsToBed(self):
        pass

    def getGTF(self):
        tmp = ''
        for feat in self.features:
            tmp += feat.toGTF()
        return tmp

    def toInterval(self):
        return intervallib.Interval(self.contig,self.start-1,self.end,self.strand)

    def fetchSequence(self,genome='hg19',connection=None):
        if connection == None:
            connection = genomelib.pygrConnect(genome)
        try:
            seq = connection[self.contig][self.start-1:self.end]
        except KeyError:
            seq = ''
        self.sequence=str(seq)
        return


#############
#lineIterator
#############
def lineIterator(gtfHandle):
    while 1:
        line = gtfHandle.readline()
        if not line: raise StopIteration
        if line.startswith("#"):continue
        gtf_entry = GTF_Entry()
        gtf_entry.read(line)
        yield gtf_entry

def GTFGeneIterator(gtfFile,verbose = False):
    handle = open(gtfFile,'r')
    iter = lineIterator(handle)
    res = {}
    if verbose:
        sys.stderr.write("Parsing GTF lines into genes...\n")
    for i in iter:
        res.setdefault(i.attributes['gene_id'],GTFGeneContainer())
        res[i.attributes['gene_id']].addFeature(i)
    for k in res.keys():
        yield res[k]

def GTFGeneIterator2(gtfFile,verbose=False):
    iter = GTFTranscriptIterator(gtfFile,verbose=verbose)
    res = {}
    for i in iter:
        res.setdefault(i.geneId,GTFGeneContainer())
        res[i.geneId].addGTFTranscript(i)
    for k in res.keys():
        yield res[k]

def GTFTranscriptIterator(gtfFile,verbose = False):
    handle = open(gtfFile,'r')
    iter = lineIterator(handle)
    res = {}
    if verbose:
        sys.stderr.write("Parsing GTF lines into transcripts...\n")
    for i in iter:
        res.setdefault(i.attributes['transcript_id'],GTFTranscriptContainer())
        res[i.attributes['transcript_id']].addFeature(i)
    for k in res.keys():
        yield res[k]

def GTFAttributeDict(gtfFile):
    """Returns a dictionary of attributes for each unique gene_id"""
    handle = open(gtfFile)
    res = {}
    fields = set([])
    for line in handle:
        if line.startswith("#"):continue
        attributes = line.rstrip().split("\t")[8].split(";")[:-1]
        attrs = [ x.strip().split(" ")[0] for x in attributes]
        fields.update(attrs)
        values = [ x.strip().split(" ")[1].strip('"') for x in attributes]
        myDict = dict(zip(attrs,values))
        res.setdefault(myDict['gene_id'],{})
        for k,v in myDict.iteritems():
            res[myDict['gene_id']].setdefault(k,set([])).add(v)
    return res

def GTFAttributeTable(gtfFile,outfile):
    """writes a table of attributes for each unique gene_id"""
    handle = open(gtfFile)
    outHandle = open(outfile,'w')
    res = {}
    fields = set([])
    for line in handle:
        if line.startswith("#"):continue
        attributes = line.rstrip().split("\t")[8].split(";")[:-1]
        attrs = [ x.strip().split(" ")[0] for x in attributes]
        fields.update(attrs)
        values = [ x.strip().split(" ")[1].strip('"') for x in attributes]
        myDict = dict(zip(attrs,values))
        res.setdefault(myDict['gene_id'],{})
        for k,v in myDict.iteritems():
            res[myDict['gene_id']].setdefault(k,set([])).add(v)
    
    #Print output to outHandle
    #Header
    print >>outHandle, "\t".join([str(x) for x in fields])
    
    for key in res.keys():
        outString = ''
        for field in fields:
            try:
                outString += ",".join(res[key][field]) + "\t" 
            except KeyError:
                outString += "-\t"
        outString.rstrip("\t")
        print >>outHandle, outString  
    return

def test():
    """
from RNASeq import newGTF
fname = 'linc_catalog.gtf'
iter = newGTF.GTFGeneIterator(fname)
for i in iter:
    print i.getGTF(),
    """
    pass

