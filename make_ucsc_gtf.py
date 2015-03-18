#!/usr/bin/env python
# encoding: utf-8
"""
make_ucsc_gtf.py

Created by Cole Trapnell on 2012-05-31.
Copyright (c) 2012 Cole Trapnell. All rights reserved.
"""

import sys
import getopt
import subprocess
import errno
import os
import tempfile
import warnings
import shutil
import copy
from datetime import datetime, date, time
import MySQLdb
import types


help_message = '''
Attaches descriptors from a UCSC table dump.

example:
./fix_ucsc.py ucsc_genes.gtf ucsc_descriptors.txt
'''

output_dir = "./"
logging_dir = output_dir + "logs/"
run_log = None
run_cmd = None

tmp_dir = output_dir + "tmp/"
bin_dir = sys.path[0] + "/"

fail_str = "\t[FAILED]\n"

# Used for prepending to cuffcompare ids to manage name collisions
ucsc_prefix_str = "ucsc_"

genome_fasta = None

class Usage(Exception):
    def __init__(self, msg):
        self.msg = msg

def right_now():
    curr_time = datetime.now()
    return curr_time.strftime("%c")

def get_version():
   return "0.0.1"

def prepare_output_dir():

    print >> sys.stderr, "[%s] Preparing output location %s" % (right_now(), output_dir)
    if os.path.exists(output_dir):
        pass
    else:        
        os.mkdir(output_dir)

    if os.path.exists(logging_dir):
        pass
    else:        
        os.mkdir(logging_dir)

    if os.path.exists(tmp_dir):
        pass
    else:        
        os.mkdir(tmp_dir)

def formatTD(td):
    hours = td.seconds // 3600
    minutes = (td.seconds % 3600) // 60
    seconds = td.seconds % 60
    return '%02d:%02d:%02d' % (hours, minutes, seconds) 

def tmp_name():
    tmp_root = output_dir + "tmp/"
    if os.path.exists(tmp_root):
        pass
    else:        
        os.mkdir(tmp_root)
    return tmp_root + os.tmpnam().split('/')[-1]

def add_attributes_with_cuffcompare(prefix, ucsc_gtf):

    print >> sys.stderr, "[%s] Augmenting GTF %s with p_id and tss_id attributes" % (right_now(), ucsc_gtf)
    cmd = ["cuffcompare"]

    cmd.extend(["-o", prefix])
    
    global genome_fasta
    if genome_fasta != None:
        cmd.extend(["-s", genome_fasta])
    
    cmd.extend(["-CG","-r", ucsc_gtf, ucsc_gtf])

    
    try:       
        print >> run_log, " ".join(cmd)
        ret = subprocess.call(cmd)
        if ret != 0:
            print >> sys.stderr, fail_str, "Error: could not execute cuffcompare"
            exit(1)
    # cuffcompare not found
    except OSError, o:
        if o.errno == errno.ENOTDIR or o.errno == errno.ENOENT:
            print >> sys.stderr, fail_str, "Error: cuffcompare not found on this system.  Did you forget to include it in your PATH?"
        exit(1)

def replace_cuffcompare_ids(cuffcompare_gtf_filename, fixed_gtf_filename):
    print >> sys.stderr, "[%s] Replacing cuffcompare transcript ids" % (right_now())
    cuffcompare_gtf = open(cuffcompare_gtf_filename)
    fixed_gtf = open(fixed_gtf_filename, "w")

    global ucsc_prefix_str
    
    for line in cuffcompare_gtf:
        line = line.strip()
        cols = line.split('\t')
        if len(cols) < 8:
            continue
        
        ucsc_gene_id = None
        tss_id = None
        p_id = None
        
        attrs = cols[8]
        attr_cols = attrs.split(';')
        for col in attr_cols:
            col = col.strip()
            col_kv_fields = col.split('"')
            if len(col_kv_fields) != 3:
                continue
            key_name = col_kv_fields[0].strip()
            value = col_kv_fields[1]
            #print >> sys.stderr, (key_name, value)
            if key_name == "oId":
                ucsc_gene_id = value
            if key_name == "tss_id":
                tss_id = value
            if key_name == "p_id":
                p_id = value
        
        fixed_attrs = []

        fixed_attrs.append('gene_id "%s";' % ucsc_gene_id)
        fixed_attrs.append('transcript_id "%s";' % ucsc_gene_id)
        if tss_id != None:
            tss_id = ucsc_prefix_str + tss_id
            fixed_attrs.append('tss_id "%s";' % tss_id)
        if p_id != None:
            p_id = ucsc_prefix_str + p_id
            fixed_attrs.append('p_id "%s";' % p_id)
        
        cols[8] = " ".join(fixed_attrs)
        print >> fixed_gtf, "\t".join(cols)

def pretty_print(f, d, level=-1, maxw=0, maxh=0, gap="", first_gap='', last_gap=''):
    # depending on the type of expression, it recurses through its elements
    # and prints with appropriate indentation

    # f   is the output file stream
    # d   is the data structure
    #
    # level is the number of allowed recursive calls, the depth at which
    #       the data structure is explored
    #       default: -1 means never stop recursing early
    # maxw  is the maximum width that will be printed from the last element
    #       of the recursion (when no further recursion is possible, or
    #       the maximal depth has been reached)
    #       default: 0 means every line will be printed in its entirety, regardless
    #                of how long it may be
    # maxh  (max height) is the maximum number of elements that will be
    #       printed from a list or a dictionary, at any level or recursion
    #       default: 0 means every list or dictionary will have all its elements
    #                printed, even if it contains thousands of elements
    #
    # gap is the gap to include before every element of a list/dic/tuple
    # first_gap is the opening gap before the opening bracket, parens or curly braces
    # first_gap is the closing gap before the closing bracket, parens or curly braces

    if level == 0:
        if type(d) != types.StringType: d = `d`

        if maxw and len(d) > maxw:
            final = ifab(maxw > 20, 10, maxw/2)
            f.write(first_gap+d[:maxw-final]+'...'+d[-final:]+' (%s chars)\n' % len(d))
        else: f.write(first_gap+d+'\n')
    elif type(d) == types.ListType:
        if not d:
            f.write(first_gap+"[]\n")
            return
        # recurse on lists
        f.write(first_gap+"[\n")
        h = 0
        for el in d:
            pretty_print(f, el, level-1, maxw, maxh, gap+'   ', gap+' ->', gap+'   ')
            if maxh:
                h = h+1
                if h >= maxh and maxh<len(d):
                    f.write(gap+' -> ... (%s in list)\n'%len(d))
                    break
        f.write(last_gap+"]\n")
    elif type(d) == types.TupleType:
        if not d:
            f.write(first_gap+"()\n")
            return
        # recurse on tuples
        f.write(first_gap+"(\n")
        h = 0
        for el in d:
            pretty_print(f, el,
                         level     = level-1,
                         maxw      = maxw,
                         maxh      = maxh,
                         gap       = gap+'   ',
                         first_gap = gap+' =>',
                         last_gap  = gap+'   ')
            if maxh:
                h = h+1
                if h >= maxh and maxh<len(d):
                    f.write(gap+' => ... (%s in tuple)\n'%len(d))
                    break
        f.write(last_gap+")\n")
    elif type(d) == types.DictType:
        if not d:
            f.write(first_gap+"{}\n")
            return
        # recurse on dictionaries
        f.write(first_gap+"{\n")
        keys = d.keys()
        keys.sort()
        key_strings = map(lambda k: ifab(type(k)==types.StringType, k, `k`), keys)
        maxlen = max(map(len, key_strings))
        h = 0
        for k,key_string in map(None, keys, key_strings):
            key_string = sfill(key_string,maxlen,'.')
            blank_string = ' '*len(key_string)
            pretty_print(f, d[k],
                         level     = level-1,
                         maxw      = maxw,
                         maxh      = maxh,
                         gap       = gap+'    %s'%blank_string,
                         first_gap = gap+'  %s: '%key_string,
                         last_gap  = gap+'    %s'%blank_string)
            if maxh:
                h = h+1
                if h >= maxh and maxh<len(keys):
                    remaining_keys = []
                    for k in keys[h:]:
                        if type(k) == types.TupleType:
                            remaining_keys.append(`k`)
                        else:
                            remaining_keys.append('%s'%k)
                    remaining_keys = string.join(remaining_keys,',')
                    #f.write(gap+'  %s (%s keys)\n'%(remaining_keys, len(keys)))
                    pretty_print(f, '  %s (%s keys)'%(remaining_keys, len(keys)),0,maxw,0,
                                 gap,gap,'')
                    break

            #gap+' '*(len(key_string)+3), '', gap+' '*(len(key_string)+5))
        f.write(last_gap+"}\n")
    elif type(d) == types.InstanceType:
        fields = dir(d)

        if not fields:
            f.write(first_gap+"*EmptyClass*\n")
            return
        # recurse on classes
        f.write(first_gap+"*ClassInstance %s\n"%d)
        fields.sort()
        key_strings = map(lambda k: ifab(type(k)==types.StringType, k, `k`), fields)
        maxlen = max(map(len, key_strings))
        h = 0
        for k,key_string in map(None, fields, key_strings):
            key_string = sfill(key_string,maxlen,'.')
            blank_string = ' '*len(key_string)
            pretty_print(f, eval('d.'+k),
                         level     = level-1,
                         maxw      = maxw,
                         maxh      = maxh,
                         gap       = gap+'    %s'%blank_string,
                         first_gap = gap+'  %s: '%key_string,
                         last_gap  = gap+'    %s'%blank_string)
            if maxh:
                h = h+1
                if h >= maxh and maxh<len(keys):
                    remaining_keys = []
                    for k in keys[h:]:
                        if type(k) == type(()):
                            remaining_keys.append(`k`)
                        else:
                            remaining_keys.append('%s'%k)
                    remaining_keys = string.join(remaining_keys,',')
                    #f.write(gap+'  %s (%s keys)\n'%(remaining_keys, len(keys)))
                    pretty_print(f,
                                 '  %s (%s keys)'%(remaining_keys, len(keys)),
                                 0,
                                 maxw,
                                 0,
                                 gap,
                                 gap,
                                 '')
                    break

            #gap+' '*(len(key_string)+3), '', gap+' '*(len(key_string)+5))
        f.write(last_gap+"*\n")
    elif type(d) == type(""):
        # simply print strings (no quotes)
        if maxw and len(d)>maxw:
            final = ifab(maxw > 20, 10, maxw/2)
            f.write(first_gap+d[:maxw-final]+'..'+d[-final:]+' (%s)\n' % len(d))
        else:
            f.write(first_gap+d+'\n')
    else:
        # string conversion of all other types
        if maxw and len(`d`)>maxw:
            final = ifab(maxw > 20, 10, maxw/2)
            f.write(first_gap+`d`[:maxw-final]+'..'+`d`[-final:]+' (%s)\n' % len(`d`))
        else:
            f.write(first_gap+`d`+'\n')

###################################
#
#Boolean Functions
#
####################################
def ifab(test, a, b):
    """x = ifab(test, a, b)
       WARNING:  Both 'a' and 'b' are evaluated
       C equivalent: x = test?a:b;
       Scheme equiv: (set x (if test a b))
       Python equiv: test and a or b
       None of the equivalents evaluates both arguments
    """
    if test: return a
    else: return b

###################################
#
#String Functions
#
####################################
def sfill(s, length, fill_char = '.'):
    #  Appends fill_char to the string s until it reaches length length
    #  ex:  sfill('hello',18,'.') -> hello...............
    #                                <---  18 chars  --->
    # useful for printing dictionaries in a cute way
    #    one......: 1
    #    five.....: 5
    #    seventeen: 17


    #list = map(None, s)
    #list.extend(map(None, fill_char*(length - len(list))))
    #return string.join(list, '')

    return s + fill_char*(length-len(s))

def rstrips(s, suffix):
    if suffix and s.endswith(suffix):
        s = s[:-len(suffix)]
    return s

def pp(d,level=-1,maxw=0,maxh=0,parsable=0):
    """ wrapper around pretty_print that prints to stdout"""
    if not parsable:
        pretty_print(sys.stdout, d, level, maxw, maxh, '', '', '')
    else:
        import pprint
        if maxw: pp2 = pprint.PrettyPrinter(width=maxw, indent=1)#, depth=level
        else: pp2 = pprint.PrettyPrinter(indent=1)#, depth=level
        pp2.pprint(d)

#################
# UCSC Genome Browser Connection options
#################

def gbdbConnect(gbdbname = "mm9"):
    gbHost = "genome-mysql.cse.ucsc.edu"
    gbUser = "genome"
    gbdb = MySQLdb.connect(host=gbHost,user=gbUser,db=gbdbname)
    return gbdb.cursor(MySQLdb.cursors.DictCursor)

def valorGbdbConnect(gbdbname='mm9'):
    gbHost = 'localhost'
    gbUser = 'root'
    gbPass = 'rinnlab1'
    gbdb = MySQLdb.connect(host=gbHost,user=gbUser,passwd=gbPass,db=gbdbname)
    return gbdb.cursor(MySQLdb.cursors.DictCursor)

########################
# Querying for IDs
########################

def fetchUCSCAnnot(genome = 'mm9',testing=False):
    """UCSC genes annotation retrieval"""
    cursor=gbdbConnect(gbdbname=genome)
    #select="SELECT * FROM refGene"
    select = """SELECT kg.name,kg.proteinID,kgti.category,kgti.genoMapCount,kgti.startComplete,kgti.endComplete,kgti.nonsenseMediatedDecay,kgx.mRNA,kgx.geneSymbol,kgx.description,ki.clusterId,ktl.value FROM knownGene kg LEFT JOIN kgTxInfo kgti ON kg.name=kgti.name LEFT JOIN kgXref kgx on kg.name=kgx.kgID LEFT JOIN knownIsoforms ki ON kg.name=ki.transcript LEFT JOIN knownToLocusLink ktl ON kg.name=ktl.name"""
    if testing:
        select += " LIMIT 10"
    cursor.execute(select)
    rows=cursor.fetchall()
    sys.stderr.write("%d rows returned\n" % len(rows))
    res = {}
    for row in rows:
        res[row['name']] = row
    sys.stderr.write("%d unique UCSC IDs\n" % len(res))
    return res

def print_UCSC_tab_delimited(UCSCDict, desc_file):
    print >> desc_file, "transcript_id\tbiotype\tgene_name\tstart_codon\tend_codon\taccession\tnum_genome_aligns\tnonsense_mediated_decay\tentrez_id\tdescription"
    for key, attrs in UCSCDict.iteritems():
        att_list = (attrs.get("category"),
                    attrs.get("geneSymbol"),
                    attrs.get("startComplete"),
                    attrs.get("endComplete"),
                    attrs.get("mRNA"),
                    attrs.get("genoMapCount"),
                    attrs.get("nonsenseMediatedDecay"),
                    attrs.get("value"),
                    attrs.get("description"))
        att_list = [str(x) if x != None else "-" for x in att_list]
        print >> desc_file, "%s\t%s" % (key, "\t".join(att_list))

def main(argv=None):
    if argv is None:
        argv = sys.argv
    try:
        try:
            opts, args = getopt.getopt(argv[1:], "ho:vs:", ["help", "output=", "genome-fasta"])
        except getopt.error, msg:
            raise Usage(msg)
    
        # option processing
        for option, value in opts:
            if option == "-v":
                verbose = True
            if option in ("-h", "--help"):
                raise Usage(help_message)
            if option in ("-s", "--genome-fasta"):
                global genome_fasta
                genome_fasta = value
            if option in ("-o", "--output-dir"):
                global output_dir
                global logging_dir
                global tmp_dir
                output_dir = value + "/"
                logging_dir = output_dir + "logs/"
                tmp_dir = output_dir + "tmp/"
        
        start_time = datetime.now()
        prepare_output_dir()
        
        global run_log
        run_log = open(logging_dir + "run.log", "w", 0)
        global run_cmd
        run_cmd = " ".join(argv)
        print >> run_log, run_cmd
                
        if len(args) < 1:
            raise Usage(help_message)
        ucsc_gtf_filename = args[0]
        ucsc_gtf = open(ucsc_gtf_filename)
        
        UCSCDict = fetchUCSCAnnot('mm9',testing=False)
        
        ucsc_desc_filename = output_dir + "/ucsc_desc.txt"
        ucsc_desc = open(ucsc_desc_filename, "w")
        
        print_UCSC_tab_delimited(UCSCDict, ucsc_desc)
        
        cuffcompare_prefix = tmp_dir + "/cuffcmp"
        add_attributes_with_cuffcompare(cuffcompare_prefix, ucsc_gtf_filename)
        
        id_replaced_filename = cuffcompare_prefix + ".id_fixed.combined.gtf"
        replace_cuffcompare_ids(cuffcompare_prefix + ".combined.gtf", id_replaced_filename)
        id_replace_file = open(id_replaced_filename)
        
        #########################
        
        #pp(UCSCDict)
        
        # print mmUCSCDict
        # 
        # header = ucsc_desc.readline()
        # 
        # id_to_rna_type = {}
        # for line in ucsc_desc:
        #     line = line.strip()
        #     cols = line.split('\t')
        #     
        #     if len(cols) < 9:
        #         continue
        #     ucsc_id = cols[0]
        #     rna_type = cols[2]    
        #     accession = cols[7]
        #     gene_symbol = cols[8]
        #     
        #     id_to_rna_type[ucsc_id] = (rna_type, accession, gene_symbol)
            
        ###########################
        
        fixed_ucsc_gtf_filename = output_dir + "/fixed_ucsc.gtf"
        fixed_ucsc_gtf = open(fixed_ucsc_gtf_filename, "w")
        
        #print id_to_rna_type
        for line in id_replace_file:
            line = line.strip()
            cols = line.split('\t')
            if len(cols) < 8:
                continue
            #print line
            rnatype = None
            accession = None    
            gene_symbol = None
            ucsc_gene_id = None
            ucsc_transcript_id = None
            ucsc_protein_id = None
            tss_id = None
            p_id = None
            
            attrs = cols[8]
            attr_cols = attrs.split(';')
            for col in attr_cols:
                col = col.strip()
                col_kv_fields = col.split('"')
                if len(col_kv_fields) != 3:
                    continue
                key_name = col_kv_fields[0].strip()
                value = col_kv_fields[1]
                #print >> sys.stderr, (key_name, value)
                if key_name == "gene_id":
                    ucsc_gene_id = value
                elif key_name == "transcript_id":
                    #metadata = id_to_rna_type.get(value)
                    metadata = UCSCDict.get(value)
                    if metadata == None:
                        print >> sys.stderr, "No metadata for %s" % value
                        continue
                    
                    rna_type = metadata.get("category")
                    accession = metadata.get("mRNA")
                    gene_symbol = metadata.get("geneSymbol")
                    ucsc_protein_id = metadata.get("proteinID")
                    ucsc_transcript_id = value
                    
                if key_name == "tss_id":
                    tss_id = value
                if key_name == "p_id":
                    p_id = value
                
            
            #print (rnatype, accession, gene_symbol, ucsc_transcript_id)        
            # Fix up the old record with the metadata we've collected        
            if rna_type != None and rna_type != "":
                cols[1] = rna_type
            
            fixed_attrs = []
            if gene_symbol != None and gene_symbol != "":
                fixed_attrs.append('gene_id "%s";' % gene_symbol)
                fixed_attrs.append('gene_name "%s";' % gene_symbol)
            else:
                fixed_attrs.append('gene_id "%s";' % ucsc_gene_id)
                
            fixed_attrs.append('transcript_id "%s";' % ucsc_transcript_id)
            
            if tss_id != None:
                fixed_attrs.append('tss_id "%s";' % tss_id)
            if ucsc_protein_id != None and ucsc_protein_id == "":
                if p_id != None:
                    fixed_attrs.append('p_id "%s";' % p_id)
            else:
                fixed_attrs.append('p_id "%s";' % ucsc_protein_id)
                
            # if accession != None:
            #                 fixed_attrs.append('transcript_id "%s";' % accession)
            #             else:
            
            
            cols[8] = " ".join(fixed_attrs)
            print >> fixed_ucsc_gtf, "\t".join(cols)
        
        finish_time = datetime.now()
        duration = finish_time - start_time
        run_log.flush()
        run_log.close()
        print >> sys.stderr,"-----------------------------------------------"
        print >> sys.stderr, "GTF construction complete [%s elapsed]" %  formatTD(duration)
                
    except Usage, err:
        print >> sys.stderr, sys.argv[0].split("/")[-1] + ": " + str(err.msg)
        print >> sys.stderr, "\t for help use --help"
        return 2


if __name__ == "__main__":
    sys.exit(main())
