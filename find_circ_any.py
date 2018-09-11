#!/usr/bin/env python
# The codes are adapted from find_circ.py to identify circRNAs from 
# any genomic loci without restriction of splicing signal.
import os
import sys
import re
import mmap
import pysam
import numpy
from logging import debug,warning,error,getLogger
from optparse import OptionParser

COMPLEMENT = {
    'a' : 't',
    't' : 'a',
    'c' : 'g',
    'g' : 'c',
    'k' : 'm',
    'm' : 'k',
    'r' : 'y',
    'y' : 'r',
    's' : 's',
    'w' : 'w',
    'b' : 'v',
    'v' : 'b',
    'h' : 'd',
    'd' : 'h',
    'n' : 'n',
    'A' : 'T',
    'T' : 'A',
    'C' : 'G',
    'G' : 'C',
    'K' : 'M',
    'M' : 'K',
    'R' : 'Y',
    'Y' : 'R',
    'S' : 'S',
    'W' : 'W',
    'B' : 'V',
    'V' : 'B',
    'H' : 'D',
    'D' : 'H',
    'N' : 'N',
}

def complement(s):
    return "".join([COMPLEMENT[x] for x in s])

def rev_comp(seq):
    return complement(seq)[::-1]

class mmap_fasta(object):
    def __init__(self,fname):
        f = file(fname)
        header = f.readline()
        row = f.readline()

        self.ofs = len(header)
        self.lline = len(row)
        self.ldata = len(row.strip())
        self.skip = self.lline-self.ldata
        self.skip_char = row[self.ldata:]
        #print "SKIP",self.skip,self.skip_char
        self.mmap = mmap.mmap(f.fileno(),0,prot=mmap.PROT_READ)

    def __getslice__(self,start,end):
        l_start = start / self.ldata
        l_end = end / self.ldata
        #print "lines",l_start,l_end
        ofs_start = l_start * self.skip + start + self.ofs
        ofs_end = l_end * self.skip + end + self.ofs
        #print "ofs",ofs_start,ofs_end
        
        s = self.mmap[ofs_start:ofs_end].replace(self.skip_char,"")
        L = end-start
        if len(s) == L:
            return s
        else:
            return s+"N"*(L-len(s))
        return 

class Accessor(object):
    supports_write = False

    def get_data(self,start,end,sense):
        return []

    def get_oriented(self,start,end,sense):
        data = self.get_data(start,end,sense)
        if sense == "-": #  self.sense_specific and
            return data[::-1]
        else:
            return data

    def get_sum(self,start,end,sense):
        return self.get_data(start,end,sense).sum()

    def flush(self):
        pass

def to_bool(obj):
    if obj == "False":
        return False
    else:
        return bool(obj)
        
class Track(object):
    """
    Abstraction of chromosome-wide adressable data like sequences, coverage, scores etc.
    Actual access to the data is delegated to accessor objects which are instantiated on-the-fly for
    each chromosome (strand) upon first access and then cached.
    Use of mmap for the accessors is recommended and implemented for sequences and numpy (C-type)
    arrays.

    See io/track_accessors.py for more examples.
    """

    def __init__(self,path,accessor,sense_specific=True,description="unlabeled track",system="hg19",dim=1,auto_flush=False,mode="r",**kwargs):
        self.path = path
        self.mode = mode
        self.acc_cache = {}
        self.accessor = accessor
        self.kwargs = kwargs
        self.sense_specific = to_bool(sense_specific)
        self.dim = int(dim)
        self.description = description
        self.auto_flush = auto_flush
        self.last_chrom = ""
        self.logger = getLogger("Track('%s')" % path)
        self.system = system

        self.logger.debug("Track(auto_flush=%s)" % (str(auto_flush)))
        kwargs['sense_specific'] = self.sense_specific
        kwargs['mode'] = self.mode
        kwargs['system'] = self.system
        kwargs['description'] = self.description
        kwargs['system'] = self.system
        kwargs['dim'] = self.dim

        if "w" in mode:
            # make sure the path exists right away so that the accessors 
            # can flush the actual data there!
            from sequence_data.io import ensure_path
            ensure_path(self.path)

    def load(self,chrom,sense):

        # automatically flush buffers whenever a new chromosome is seen. reduces memory-footprint for sorted input data
        if self.auto_flush and chrom != self.last_chrom:
            self.logger.debug("Seen new chromosome %s. Flushing accessor caches." % chrom)
            self.flush_all()
        self.last_chrom = chrom

        ID = self.get_identifier(chrom,sense)
        if not ID in self.acc_cache:
            self.logger.debug("Cache miss for %s%s. creating new accessor" % (chrom,sense))
            self.acc_cache[ID] = self.accessor(self.path,chrom,sense,**(self.kwargs))

        return self.acc_cache[ID]

    def save(self):
        """
        If a track is opened in "rw" or "w" mode this will save the track-definition config files and flush all accessors.
        Saving of the actual data is performed by the accessors that support writing upon a call to flush.
        """
        if not "w" in self.mode:
            self.logger.warning("save() called on a read-only opened track. Ignored!")
            return

        if not self.accessor.supports_write:
            self.logger.warning("save() called on a track with only read-access supporting accessors. Ignored!")
            return
      
        self.logger.debug("save(): writing '%s'" % self.path)

        def to_str(obj):
            # convert simple data-types to their string representation
            # but classes and more complex types to their names.
            return getattr(obj,"__name__",str(obj))

        kwarg_str = "\n".join(["%s=%s" % (k,to_str(self.kwargs[k])) for k in sorted(self.kwargs.keys()) if k != "mode"])
        file(os.path.join(self.path,"track.rc"),"w+").write(trackrc % dict(accessor=self.accessor.__name__,kwargs=kwarg_str))
        self.flush_all()

    def __del__(self):
        if "w" in self.mode:
            self.save()

    def flush(self,chrom,sense):
        ID = self.get_identifier(chrom,sense)
        if ID in self.acc_cache:
            self.logger.warning("Flushing %s%s" % (chrom,sense))
            del self.acc_cache[ID]

    def flush_all(self):
        for a in self.acc_cache.values():
            a.flush()
        self.acc_cache = {}

    def get(self,chrom,start,end,sense):
        acc = self.load(chrom,sense)
        return acc.get_data(start,end,sense)

    def get_oriented(self,chrom,start,end,sense):
        acc = self.load(chrom,sense)
        return acc.get_oriented(start,end,sense)

    def get_sum(self,chrom,start,end,sense):
        acc = self.load(chrom,sense)
        return acc.get_sum(start,end,sense)
        
    def get_identifier(self,chrom,sense):
        if self.sense_specific:
            return chrom+sense
        else:
            return chrom

class GenomeAccessor(Accessor):
    def __init__(self,path,chrom,sense,**kwargs):
        debug("# GenomeAccessor mmap: Loading genomic sequence for chromosome %s from '%s'" % (chrom,path))

        fname = os.path.join(path,chrom+".fa")
        try:
            self.data = mmap_fasta(fname)
        except IOError:
            warning("Could not access '%s'. Switching to dummy mode (only Ns)" % fname)
            self.get_data = self.get_dummy
            self.get_oriented = self.get_dummy

        # TODO: maybe remove this if not needed
        self.get = self.get_oriented

    def get_data(self,start,end,sense):
        if start < 0 or end < 0:
            return self.get_dummy(start,end,sense)
        #UCSC convention: start with 1, end is inclusive
        if sense == "+":
            return self.data[start:end]
        else:
            return complement(self.data[start:end])

    def get_oriented(self,start,end,sense):
        if end < 0:
            return self.get_dummy(start,end,sense)
        elif start < 0:
            return self.get_dummy(start,0,sense) + self.get_oriented(0,end,sense)
        if sense == "+":
            return self.data[start:end]
        else:
            try:
                return rev_comp(self.data[start:end])
            except KeyError:
                print start,end,sense

    def get_dummy(self,start,end,sense):
        return "N"*(end-start)



usage = """

  bowtie2 anchors.qfa.gz | %prog > candidates.bed 2> candidates.reads

"""

parser = OptionParser(usage=usage)
parser.add_option("-S","--system",dest="system",type=str,default="",help="model system database (alternatively use -G <path_to_genome_folder>)")
parser.add_option("-G","--genome",dest="genome",type=str,default="",help="path to genome. Point to folder with one fasta file for each chromosome.")
parser.add_option("-a","--anchor",dest="asize",type=int,default=20,help="anchor size (default=20)")
parser.add_option("-w","--wiggle",dest="wiggle",type=int,default=2,help="maximum nts a linear splice site may be away from circ to be counted as competitor (default=2)")
parser.add_option("-m","--margin",dest="margin",type=int,default=2,help="maximum nts the BP is allowed to reside inside an anchor (default=2)")
parser.add_option("-d","--maxdist",dest="maxdist",type=int,default=2,help="maximum mismatches (no indels) allowed in anchor extensions (default=2)")
parser.add_option("-p","--prefix",dest="prefix",default="",help="prefix to prepend to each name")
parser.add_option("-q","--min_uniq_qual",dest="min_uniq_qual",type=int,default=2,help="minimal uniqness (alignment score margin to next-best hit) for anchor alignments (default=2)")
parser.add_option("-s","--stats",dest="stats",default="runstats.log",help="where to put stats (default='runstats.log'")
parser.add_option("-r","--reads2samples",dest="reads2samples",default="",help="path to tab-separated two-column file with read-name prefix -> sample ID mapping")

options,args = parser.parse_args()

if options.system:
    # you have the rest of the library installed and set up, great
    import sequence_data.systems
    genome = getattr(sequence_data.systems,options.system).genome
else:
    genome = Track(options.genome,GenomeAccessor)


from itertools import izip_longest
def grouper(n, iterable, fillvalue=None):
    "grouper(3, 'ABCDEFG', 'x') --> ABC DEF Gxx"
    args = [iter(iterable)] * n
    return izip_longest(fillvalue=fillvalue, *args)

if options.reads2samples:
    samples = [line.rstrip().split('\t') for line in file(options.reads2samples)]
else:
    samples = []
    
samples.append(('','unknown'))

        
minmapscore = options.asize * (-2)
from collections import defaultdict

class Hit(object):
    def __init__(self):
        self.reads = []
        self.uniq = set()
        self.mapquals = []
        self.tissues = defaultdict(int)
        self.edits = []
        self.overlaps = []
        self.n_hits = []
        
    def add(self,read,A,B,dist,ov,n_hits):
        self.reads.append(read)
        self.edits.append(dist)
        self.overlaps.append(ov)
        self.n_hits.append(n_hits)

        # Alignment Score - Secondbest hit score
        aopt = dict(A.tags)
        bopt = dict(B.tags)
        qA = aopt.get('AS') - aopt.get('XS',minmapscore)
        qB = bopt.get('AS') - bopt.get('XS',minmapscore)
       
        self.mapquals.append((qA+qB,qA,qB))
        
        for (prefix,tiss) in samples:
            if A.qname.startswith(prefix):
                self.tissues[tiss] += 1
                break
    
        self.uniq.add((read,tiss))
        self.uniq.add((rev_comp(read),tiss))


    def scores(self,chrom,start,end,sense):
        n_reads = len(self.reads)
        n_uniq = len(self.uniq) / 2
        
        total_mq,best_qual_A,best_qual_B = sorted(self.mapquals,reverse=True)[0]

        wiggle = numpy.arange(-options.wiggle,options.wiggle+1)   
        
        spliced_at_begin = 0
        for x in wiggle:
            begin = (chrom,start+x,sense)
            if begin in loci:
                for l in loci[begin]:
                    spliced_at_begin += len(l.reads)

        spliced_at_end = 0
        for x in wiggle:
            ending = (chrom,end+x,sense)
            if ending in loci:
                for l in loci[ending]:
                    spliced_at_end += len(l.reads)

        #print self.edits,self.overlaps,self.n_hits
        tissues = sorted(self.tissues.keys())
        tiss_counts = [str(self.tissues[k]) for k in tissues]
        return (n_reads,n_uniq,best_qual_A,best_qual_B,spliced_at_begin,spliced_at_end,tissues,tiss_counts,min(self.edits),min(self.overlaps),min(self.n_hits))
                
        
loci = defaultdict(list)
circs = defaultdict(Hit)
splices = defaultdict(Hit)

N = defaultdict(int)

from numpy import chararray as carray
from numpy import fromstring,byte

def find_breakpoints(A,B,read,chrom,margin=options.margin,maxdist=options.maxdist):

    def mismatches(a,b):
        a,b = fromstring(a,dtype=byte), fromstring(b,dtype=byte)
        return (a != b).sum()
        
    L = len(read)
    hits = []
    eff_a = options.asize-margin
    internal = read[eff_a:-eff_a].upper()
        
    flank = L - 2*eff_a + 2

    A_flank = genome.get(chrom,A.aend-margin,A.aend-margin + flank,'+').upper()
    B_flank = genome.get(chrom,B.pos - flank+margin,B.pos+margin,'+').upper()
    l = L - 2*eff_a
    for x in range(l+1):
        spliced = A_flank[:x] + B_flank[x+2:]
        dist = mismatches(spliced,internal)        

        ov = 0
        if x < margin:
            ov = margin-x
        if l-x < margin:
            ov = margin-(l-x)
        
        if dist <= maxdist:
            gt = A_flank[x:x+2]
            ag = B_flank[x:x+2]
            
            start,end = B.pos+margin-l+x,A.aend-margin+x+1
            start,end = min(start,end),max(start,end)
            if gt == 'GT' and ag == 'AG':
                # GT/AG is used as a character not a constraint
                splice = 0
                hits.append((dist,ov,splice,chrom,start,end,'+'))
            elif gt == 'CT' and ag == 'AC':
                splice = 0
                hits.append((dist,ov,splice,chrom,start,end,'-'))
            elif gt == 'AT' and ag == 'AC':
                # minor splicing signal
                splice = 1
                hits.append((dist,ov,splice,chrom,start,end,'+'))
            elif gt == 'GT' and ag == 'AT':
                splice = 1
                hits.append((dist,ov,splice,chrom,start,end,'-'))
            else:
                splice = 2
                # + strand is meaningless here
                hits.append((dist,ov,splice,chrom,start,end,'+'))


    if len(hits) < 2:
        return hits

    # return only hits that are tied with the best candidate by edit distance and anchor overlap. 
    # Hits are still sorted, with low edit distance beating low anchor overlap
    hits = sorted(hits)
    best = hits[0]
    return [h for h in hits if ((h[0] == best[0]) and (h[1] == best[1]))]
    
    
sam = pysam.Samfile(sys.argv[1],'r')
try:
    for A,B in grouper(2,sam):

        N['total'] += 1
        if A.is_unmapped or B.is_unmapped:
            N['unmapped'] += 1
            continue
        if A.tid != B.tid:
            N['other_chrom'] += 1
            continue
        if A.is_reverse != B.is_reverse:
            N['other_strand'] += 1
            continue

        dist = B.pos - A.pos
        if numpy.abs(dist) < options.asize:
            N['overlapping_anchors'] += 1
            continue
        
        if (A.qname)and((A.is_reverse and dist > 0) or (not A.is_reverse and dist < 0)):
            # get original read sequence
            read0 = A.qname.split('__')[1]
            chrom = sam.getrname(A.tid)
            
            if A.is_reverse:
                #print "ISREVERSE"
                A,B = B,A
                read = rev_comp(read0)
                sense1 = '-'
                bp = find_breakpoints(A,B,read,chrom)
            else:
                bp = find_breakpoints(A,B,read0,chrom)
                sense1 = '+'

            if not bp:
                N['circ_no_bp'] += 1
            else:
                N['circ_reads'] += 1

            n_hits = len(bp) 
            for h in bp:
                #print h
                # for some weird reason for circ we need a correction here
                dist,ov,splice,chrom,start,end,sense2 = h
                if splice == 0 and sense1 == sense2:
                    splice = 2
                elif splice == 1 and sense1 == sense2:
                    splice = 1
                else:
                    splice = 0
                h = (chrom,start+1,end-1,sense1,splice)
                circs[h].add(read0,A,B,dist,ov,n_hits)

        if (A.qname)and((A.is_reverse and dist < 0) or (not A.is_reverse and dist > 0)):
            read0 = A.qname.split('__')[1]
            chrom = sam.getrname(A.tid)
            
            if A.is_reverse:
                #print "ISREVERSE"
                A,B = B,A
                read = rev_comp(read0)
                sense1 = '-'
                bp = find_breakpoints(A,B,read,chrom)
            else:
                bp = find_breakpoints(A,B,read0,chrom)
                sense1 = '+'
            if not bp:
                N['splice_no_bp'] += 1
            else:
                N['spliced_reads'] += 1
            n_hits = len(bp)                
            for h in bp:
                #print h
                dist,ov,splice,chrom,start,end,sense2 = h
                if splice == 0 and sense1 == sense2:
                    splice = 2
                elif splice == 1 and sense1 == sense2:
                    splice = 1
                else:
                    splice = 0
                h = (chrom,start,end,sense1,splice)
                splices[h].add(read0,A,B,dist,ov,n_hits)
                loci[(chrom,start,sense1)].append(splices[h])
                loci[(chrom,end,sense1)].append(splices[h])

            
    #break
except KeyboardInterrupt:
    pass

def output(ofile,rfile,cand,prefix):
    head="\t".join(['chrom','start','end','name','n_reads','strand','n_uniq','splice','best_qual_A','best_qual_B',
        'spliced_at_begin','spliced_at_end','tissues','tiss_counts','edits','anchor_overlap','breakpoints'])
    ofile.write("#%s\n" % head)

    #print "#","\t".join(['chrom','start','end','name','n_reads','strand','n_uniq','best_qual_A','best_qual_B','spliced_at_begin','spliced_at_end','tissues','tiss_counts','edits','anchor_overlap','breakpoints'])
    n = 1
    for c,hit in cand.items():
        #print c
        chrom,start,end,sense,splice = c
        n_reads,n_uniq,best_qual,best_other,spliced_at_begin,spliced_at_end,tissues,tiss_counts,min_edit,min_anchor_ov,n_hits = hit.scores(chrom,start,end,sense)
        
        if best_other < options.min_uniq_qual:
            N['anchor_not_uniq'] += 1
            continue
        
        name = "%s_%s_%06d" % (options.prefix,prefix,n)
        n += 1
        #sys.stderr.write("%s\t%s\n" % (name,"\t".join(sorted(reads))))
        for r in hit.reads:
            rfile.write("%s\t%s\n" % (name,r))
        
        bed = [chrom,start-1,end,name,n_reads,sense,n_uniq,splice,best_qual,best_other,spliced_at_begin,spliced_at_end,",".join(tissues),",".join(tiss_counts),min_edit,min_anchor_ov,n_hits]
        item="\t".join([str(b) for b in bed])
        ofile.write("%s\n" % item)
        #print "\t".join([str(b) ffor b in bed])

stats = file(options.stats,"w")
stats.write(str(N)+"\n")
ofile = open(sys.argv[2],'w')
rfile = open(sys.argv[3],'w')

output(ofile,rfile,circs,"circ")
output(ofile,rfile,splices,"norm")

