#!/usr/bin/python
import sys,types,string
#############
#pygr tools
#############
class Annot:
    """Annotation class for pygr data"""
    def __init__(self,name,chr,strand,start,end):
        self.name=name
        self.chr=chr
        self.strand=strand
        self.start=start
        self.end=end
        
##################
#nuID implementation for python
###################
def mreplace(s,chararray=['A','C','G','T','U'],newarray=['0','1','2','3','3']):
    for a,b in zip(chararray,newarray):
        s=s.replace(a,b)
    return s

def seq2nuID(seq):
    """Converts a string DNA or RNA sequence into its corresponding 'nuID'"""
    
    """ 
        Default code includes "_" as char.  This conflicts with parsing for shrimp.  So for my specific instance,
        "_" has been replaced with "!"
    """
    code = map(chr,range(65,91))+map(chr,range(97,123))+map(str,range(0,10))+map(str,("!","."))
    seq=seq.upper()
    num=mreplace(seq)
    if len(num)%3!=0:
        appLen = 3-len(num)%3
        num = num+(appLen*str(0))
    else:
        appLen=0
    numArray=[]
    checkSumArray=[]
    charArray=[]
    for i in range(0,len(num),3):
        subnum=num[i:i+3]
        numArray.append(subnum)
        code64=int(subnum[0])*4**2+int(subnum[1])*4+int(subnum[2])
        checkSumArray.append(code64)
        charArray.append(code[code64])
    checkSum=sum(checkSumArray)
    res=checkSum%21
    checkCode=code[res*3+appLen]
    #print numarray
    #print checkSumArray
    id = str(checkCode)+"".join(charArray)
    return id

def nuID2seq(nuID):
    """ 
        Default code includes "_" as char.  This conflicts with parsing for shrimp.  So for my specific instance,
        "_" has been replaced with "!"
    """
    import math
    code = map(chr,range(65,91))+map(chr,range(97,123))+map(str,range(0,10))+map(str,("!","."))
    ind=range(1,len(code)+1)
    names=dict(zip(code,ind))
    numArray=[]
    for l in nuID:
        numArray.append(names[l]-1) #Possibly need to add -1?
    checkCode=int(numArray.pop(0))
    if checkCode==63:
        assert "Coding error or not a nuID!\nCheck code should not include '.'!"
    cutlen=checkCode%3
    res=int(math.floor(checkCode/3))
    num=''
    for i in numArray:
        subDecode = [int(math.floor(i/4**2)),int(math.floor((i%4**2)/4)),int(i%4)]
        newsub = "".join(str(j) for j in subDecode)
        num=num+newsub
    checkSum=sum(numArray)
    if res != checkSum%21:
        assert "Coding Error or not a nuID"
    nucleotide=["A","C","G","T"]
    seq=mreplace(num,['0','1','2','3'],nucleotide)
    seq=seq[:-1]
    return seq

######
#
#Dictionary tools
#
#######

def sort_by_value(d):
    """ Returns the keys of dictionary d sorted by their values """
    items=d.items()
    backitems=[ [v[1],v[0]] for v in items]
    backitems.sort(reverse=True)
    return [ backitems[i][1] for i in range(0,len(backitems))]

def sbv2(d,reverse=False):  
    ''' proposed in PEP 265, using  the itemgetter '''  
    from operator import itemgetter
    return sorted(d.iteritems(), key=itemgetter(1), reverse=True)  

def sortListofDicts(fieldname):
    """useful for sorting a list of dictionaries by a given key (fieldname)
    usage:
    mylist.sort(sortListofDicts('start')  #will sort a list of intervals by i['start']
    """
    def compare_two_dicts (a,b):
        return cmp(a[fieldname],b[fieldname])
    return compare_two_dicts

def sort_dict(d,reverse=True):
    return sorted(d.iteritems(), key=lambda (k,v): (v,k), reverse=reverse)

########
#
#Pretty Printing
#
########
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

def pp(d,level=-1,maxw=0,maxh=0,parsable=0):
    """ wrapper around pretty_print that prints to stdout"""
    if not parsable: 
        pretty_print(sys.stdout, d, level, maxw, maxh, '', '', '')
    else:
        import pprint
        if maxw: pp2 = pprint.PrettyPrinter(width=maxw, indent=1)#, depth=level
        else: pp2 = pprint.PrettyPrinter(indent=1)#, depth=level
        pp2.pprint(d)

def test_pp():
    pp({'one': ('two',3,[4,5,6]),
        7: (lambda x: 8*9),
        'ten': ['ele', {'ven': 12,
                        (13,14): '15'}]})

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

def hamming_distance(s1, s2):
    """Returns the hamming (or edit) distance between two strings or list of iterable elements"""
    assert len(s1) == len(s2)
    return sum(ch1 != ch2 for ch1, ch2 in zip(s1, s2))

######################################
#
#Ranking and Ordering
#
######################################
from random import uniform, sample

def order(x, NoneIsLast = True, decreasing = False):
    """
    Returns the ordering of the elements of x. The list
    [ x[j] for j in order(x) ] is a sorted version of x.

    Missing values in x are indicated by None. If NoneIsLast is true,
    then missing values are ordered to be at the end.
    Otherwise, they are ordered at the beginning.
    """
    omitNone = False
    if NoneIsLast == None:
        NoneIsLast = True
        omitNone = True
        
    n  = len(x)
    ix = range(n)
    if None not in x:
        ix.sort(reverse = decreasing, key = lambda j : x[j])
    else:
        # Handle None values properly.
        def key(i, x = x):
            elem = x[i]
            # Valid values are True or False only.
            if decreasing == NoneIsLast:
                return not(elem is None), elem
            else:
                return elem is None, elem
        ix = range(n)
        ix.sort(key=key, reverse=decreasing)
            
    if omitNone:
        n = len(x)
        for i in range(n-1, -1, -1):
            if x[ix[i]] == None:
                n -= 1
        return ix[:n]
    return ix


def rank(x, NoneIsLast=True, decreasing = False, ties = "first"):
    """
    Returns the ranking of the elements of x. The position of the first
    element in the original vector is rank[0] in the sorted vector.

    Missing values are indicated by None.  Calls the order() function.
    Ties are NOT averaged by default. Choices are:
         "first" "average" "min" "max" "random" "average"
    """
    omitNone = False
    if NoneIsLast == None:
        NoneIsLast = True
        omitNone = True
    O = order(x, NoneIsLast = NoneIsLast, decreasing = decreasing)
    R = O[:]
    n = len(O)
    for i in range(n):
        R[O[i]] = i
    if ties == "first" or ties not in ["first", "average", "min", "max", "random"]:
        return R
        
    blocks     = []
    isnewblock = True
    newblock   = []
    for i in range(1,n) :
        if x[O[i]] == x[O[i-1]]:
            if i-1 not in newblock:
                newblock.append(i-1)
            newblock.append(i)
        else:
            if len(newblock) > 0:
                blocks.append(newblock)
                newblock = []
    if len(newblock) > 0:
        blocks.append(newblock)

    for i, block  in enumerate(blocks):
        # Don't process blocks of None values.
        if x[O[block[0]]] == None:
            continue
        if ties == "average":
            s = 0.0
            for j in block:
                s += j
            s /= float(len(block))
            for j in block:
                R[O[j]] = s                
        elif ties == "min":
            s = min(block)
            for j in block:
                R[O[j]] = s                
        elif ties == "max":
            s =max(block)
            for j in block:
                R[O[j]] = s                
        elif ties == "random":
            s = sample([O[i] for i in block], len(block))
            for i,j in enumerate(block):
                R[O[j]] = s[i]
        else:
            for i,j in enumerate(block):
                R[O[j]] = j
    if omitNone:
        R = [ R[j] for j in range(n) if x[j] != None]
    return R

def uniqify(seq): 
    # Not order preserving 
    keys = {} 
    for e in seq: 
        keys[e] = 1 
    return keys.keys()