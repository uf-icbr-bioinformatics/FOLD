#!/usr/bin/env python

import sys
import csv
import os.path
import random
from numpy import var
from datetime import datetime
import subprocess as sp
import importlib

try:
    primer3 = importlib.import_module('primer3')
    HAS_PRIMER3 = True
except ImportError:
    HAS_PRIMER3 = False

# Constants

U1 = "ATCACCAACTACCCACACACACC"
U2 = "CTACCCCACCTTCCTCATTCTCT"

## Utilities

# Check that bedtools is available

def checkBedtools():
    try:
        return sp.check_output("bedtools --version", shell=True, stderr=sp.STDOUT).strip()
    except sp.CalledProcessError:
        return False

# MT calculation

def mt(gc, length):
    """Returns the MT of a sequence of length `len' containing `gc' C or G nucleotides.
Computation performed according to Varley et al., see: http://basic.northwestern.edu/biotools/OligoCalc.html."""
    p = -21.597097928022087
    if length < 14:
        return 2*(length - gc) + 4*gc
    else:
        return 100.5 + (41.0 * gc / length) - (820.0 / length) + p

# Reverse complement

BASES = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N',
         'a': 't', 'c': 'g', 'g': 'c', 't': 'a', 'n': 'n'}

def revcomp(s):
    bases = [ BASES[b] for b in s[::-1] ]
    return "".join(bases)

# Number of GC bases in a sequence

def ngc(seq):
    return seq.count("C") + seq.count("G") + seq.count("c") + seq.count("g")

def nlower(seq):
    nl = 0
    for c in seq:
        if c.islower():
            nl += 1
    return nl

# Find candidate positions to test

def findCandidatePositions(positions, howmany, start, direction):
    """Return `howmany' positions from the list of positions starting at `start' in direction `direction'
(1 = left to right, -1 = right to left)."""
    if not positions:
        return []

    # Find first position after `start'
    for i in range(len(positions)):
        if positions[i] >= start:
            break
    if direction == 1:
        return positions[i:i+howmany]
    else:
        p = max(i-howmany, 0)
        return positions[p:i]

# Oligo

class Oligo(object):
    start = 0
    end = 0
    length = 0
    gc = 0
    gcperc = 0.0
    mt = 0
    sequence = ""
    flags = ""

    def __init__(self, seq, start, direction, bestlen, bestmt):
        #print (start, direction, bestlen, bestmt)
        if direction == 1:
            self.start = start
            self.end = start + bestlen - 1
        else:
            self.start = start - bestlen + 1
            self.end = start
        self.length = bestlen
        self.mt = bestmt
        self.sequence = seq[self.start:self.end+1]
        self.gc = ngc(self.sequence)
        self.gcperc = 1.0 * self.gc / self.length

# Triple of oligos

class Triple(object):
    oligo1 = None
    oligo2 = None
    oligo3 = None
    size = 0
    tmspread = 100
    score = 1000

    def __init__(self, o1, o2, o3):
        self.oligo1 = o1
        self.oligo2 = o2
        self.oligo3 = o3
        self.size = o2.start - o1.start + 1
        self.tmspread = max(o1.mt, o2.mt, o3.mt) - min(o1.mt, o2.mt, o3.mt)

# Extracting sequences for regions

class Sequence(object):
    name = ""
    seq = ""
    chrom = ""
    start = 0
    end = 0
    strand = ""
    valid = True                # Set to False if no valid oligos found
    triples = []
    rndtriple = None            # One triple picked at random
    best = None                 # Best triple found so far
    bestGC = 0                  # GC content corresponding to the best triple

    def __init__(self, name, seq):
        self.name = name
        self.seq = seq
        self.triples = []

    def coordinates(self):
        return "{}:{}-{}".format(self.chrom, self.start, self.end)

    def write(self, filename, label=None, start=None, end=None):
        with open(filename, "w") as out:
            self.writes(out, label=label, start=start, end=end)

    def writes(self, out, label=None, start=None, end=None):
        out.write(">{}\n{}\n".format(label or self.name, self.seq[start:end]))

    def findVT(self):
        positions = []
        for i in range(1, len(self.seq)):
            if self.seq[i] in 'Tt' and self.seq[i-1] not in 'Tt':
                positions.append(i)
        return positions

    def randomTriple(self):
        self.rndtriple = random.choice(self.triples)
        return self.rndtriple

    def bestGC(self):
        start = self.best.oligo1.start+1
        end = self.best.oligo2.start+1
        gc = ngc(self.seq[start:end])
        self.bestGC = 1.0 * gc / (end-start)
        return self.bestGC

class SequenceManager(object):
    reference = None

    def __init__(self, reference):
        self.reference = reference

    def getSequences(self, intervals):
        seqs = {}
        proc = sp.Popen("bedtools getfasta -fi {} -bed /dev/stdin".format(self.reference), shell=True, stdin=sp.PIPE, stdout=sp.PIPE)
        for intv in intervals:
            #sys.stderr.write("{}\n".format(intv))
            proc.stdin.write("{}\t{}\t{}\n".format(intv[0], intv[1], intv[2]).encode('utf-8'))
            seqname = ">{}:{}-{}".format(intv[0], intv[1], intv[2])
            s = Sequence(seqname, None)
            s.chrom = intv[0]
            s.start = intv[1]
            s.end   = intv[2]
            s.strand = intv[3]
            seqs[seqname] = s
        proc.stdin.close()

        while True:
            name = proc.stdout.readline()
            if not name:
                break
            name = name.decode('utf-8').rstrip("\n")
            seq = proc.stdout.readline().decode('utf-8').rstrip("\n")
            s = seqs[name]
            s.seq = seq

        #sys.stderr.write("{}\n".format(seqs))
        return list(seqs.values())

# Intervals

class Interval(object):
    name = ""                   # Source of this interval, e.g. a gene name
    label = ""                  # chr:start-end
    chrom = ""
    start = 0
    end = 0

    def __init__(self, name, chrom, start, end):
        self.name = name
        self.chrom = chrom
        self.start = start 
        self.end = end
        self.label = "{}:{}-{}".format(chrom, start, end)

# Simple Gene DB

class Gene(object):
    accession = ""
    name = ""
    transcripts = []

    def __init__(self, accession, name):
        self.accession = accession
        self.name = name
        self.transcripts = []

    def addTranscript(self, accession, chrom, start, end, strand):
        """Add a new transcript with the given properties to this gene, unless 
one with the same TSS already exists."""
        for tr in self.transcripts:
            if (strand == '+' and start == tr[2]) or (strand == '-' and end == tr[3]):
                return 0          # already present, ignore
        self.transcripts.append([self.name + ":" + accession, chrom, start, end, strand])
        return 1

# Simple GTF parser

class GTFParser(object):
    filename = ""
    genes = {}
    fixChroms = False

    def __init__(self, filename):
        self.filename = filename
        self.genes = {}

    def ngenes(self):
        return len(self.genes)

    def fixChrom(self, c):
        if c.startswith("chr"):
            return c
        elif self.fixChroms:
            return "chr" + c
        else:
            return c

    def getName(self, annots):
        parts = annots.split(";")
        for p in parts:
            p = p.lstrip()
            if p.startswith("gene_name"):
                return p[10:].strip('"')
        return None

    def getNames(self, annots):
        wanted = ["gene_name", "gene_id", "transcript_id"]
        d = {}
        parts = annots.split(";")
        for p in parts:
            p = p.lstrip()
            for w in wanted:
                if p.startswith(w):
                    d[w] = p[len(w)+1:].strip('"')
        return d

    def load(self):
        bytx = {}
        sys.stderr.write("Loading transcripts from GTF file {}: ".format(self.filename))
        sys.stderr.flush()
        n = 0
        with open(self.filename, "r") as f:
            c = csv.reader(f, delimiter='\t')
            for line in c:
                if line[0][0] == '#':
                    continue
                if line[2] == "gene":
                    names = self.getNames(line[8])
                    g = Gene(names["gene_id"], names["gene_name"])
                    self.genes[g.name] = g
                    bytx[g.accession] = g
                elif line[2] == "transcript":
                    names = self.getNames(line[8])
                    gacc = names["gene_id"]
                    if gacc in bytx:
                        g = bytx[gacc]
                        n += g.addTranscript(names["transcript_id"],
                                             self.fixChrom(line[0]),
                                             int(line[3]), int(line[4]), line[6])
                        if n % 5000 == 0:
                            sys.stderr.write("\x1b[GLoading transcripts from GTF file {}: {}".format(self.filename, n))
                            sys.stderr.flush()
        sys.stderr.write("\x1b[GLoading transcripts from GTF file {}: {} genes, {} transcripts loaded.\n".format(self.filename, len(self.genes), n))

    def get(self, name):
        return self.genes[name] if name in self.genes else None

# And we get to also parse GFF almost for free...

class GFFParser(GTFParser):

    def getName(self, annots):
        parts = annots.split(";")
        for p in parts:
            p = p.lstrip()
            if p.startswith("Name="):
                return p[5:]
        return None

    def getNames(self, annots):
        wanted = ["Name=", "ID=gene:", "ID=transcript:", "Parent=gene:"]
        d = {}
        parts = annots.split(";")
        for p in parts:
            p = p.lstrip()
            for w in wanted:
                if p.startswith(w):
                    d[w] = p[len(w):]
        if "Name=" in d:
            d["gene_name"] = d["Name="]
        if "ID=gene:" in d:
            d["gene_id"] = d["ID=gene:"]
        if "ID=transcript:" in d:
            d["transcript_id"] = d["ID=transcript:"]
        if "Parent=gene:" in d:
            d["gene_id"] = d["Parent=gene:"]
        return d

# And with a small change we also parse refFlat files.

class refFlatParser(GTFParser):

    def load(self):
        sys.stderr.write("Loading transcripts from refFlat file {}: ".format(self.filename))
        sys.stderr.flush()
        n = 0
        with open(self.filename, "r") as f:
            c = csv.reader(f, delimiter='\t')
            for line in c:
                if line[0][0] == '#':
                    continue
                name = line[0]
                txacc = line[1]
                chrom = self.fixChrom(line[2])
                if "_" in chrom:
                    continue
                strand = line[3]
                start = int(line[4])
                end = int(line[5])
                if name in self.genes:
                    g = self.genes[name]
                else:
                    g = Gene(name, name)
                    self.genes[name] = g
                n += g.addTranscript(txacc, chrom, start, end, strand)
                if n % 5000 == 0:
                    sys.stderr.write("\x1b[GLoading transcripts from refFlat file {}: {}".format(self.filename, n))
                    sys.stderr.flush()
        sys.stderr.write("\x1b[GLoading transcripts from refFlat file {}: {} genes, {} transcripts loaded.\n".format(self.filename, len(self.genes), n))

class BEDfileParser(GTFParser):

    def load(self):
        sys.stderr.write("Loading regions from BED file {}: ".format(self.filename))
        sys.stderr.flush()
        n = 0
        with open(self.filename, "r") as f:
            c = csv.reader(f, delimiter='\t')
            for line in c:
                if not line:
                    continue
                if not(line[0]) or line[0][0] == '#':
                    continue
                name = line[3]
                txacc = name
                chrom = self.fixChrom(line[0])
                start = int(line[1])
                end = int(line[2])
                if name in self.genes:
                    g = self.genes[name]
                else:
                    g = Gene(name, name)
                    self.genes[name] = g
                n += g.addTranscript(txacc, chrom, start, end, "+")
        sys.stderr.write("\x1b[GLoading regions from BED file {}: {} regions loaded.\n".format(self.filename, len(self.genes), n))

## Reporter

class Report(object):
    out = None

    def __init__(self, out):
        self.out = out
        self.out.write("### FenGC Report - started {} ###\n\n".format(datetime.now()))

    def badGenes(self, badgenes):
        self.out.write("*** The following {} genes were not found in the gene database:\n\n".format(len(badgenes)))
        for b in badgenes:
            self.out.write(b + "\n")
        self.out.write("===\n\n")

    def badOligos(self, badoligos):
        self.out.write("*** No valid oligos could be found for the following {} transcripts:\n\n".format(len(badoligos)))
        for b in badoligos:
            self.out.write(b + "\n")
        self.out.write("===\n\n")
        
    def heterodimers(self, hets):
        sizes = [5, 6, 5, 6]
        for h in hets:
            for i in range(4):
                sizes[i] = max(sizes[i], len(h[i]))
        fmt = "{{:{}}}  {{:{}}}  {{:{}}}  {{:{}}}  ".format(*sizes)
        fmt1 = fmt + "Tm\n"
        fmt2 = fmt + "{:.2f}\n"
        self.out.write("*** {} potential heterodimers found:\n\n".format(len(hets)))
        self.out.write(fmt1.format("Gene1", "Oligo1", "Gene2", "Oligo2"))
        for h in hets:
            self.out.write(fmt2.format(*h))
        self.out.write("===\n\n")

## Main

class Main(object):
    reference = None            # Name of reference file
    genelist = None             # -g Name of genes DB
    genes = []
    sequences = []              # List of Sequence objects
    outfile = "/dev/stdout"     # -o
    reportfile = "report.txt"   # -r
    seqman = None
    toFasta = False             # -F, -f
    oligoFasta = False          # -O
    local = False               # -l If true, skips global optimization.
    transcripts = "all"         # "all" or "longest" - set by -L
    bedmode = False             # True when input is from a BED file (regions)
    fixChroms = False           # -c  If True, add "chr" in front of chromosome names in annotations file, if missing

    upstream = 400              # -u
    downstream = 100            # -d
    field = 2000                # -s
    includeTSS = True           # -nt to disable
    minlength = 8               # -ol
    maxlength = 30              # -ol
    minmt = 62                  # -tl
    maxmt = 68                  # -th
    mingc = 40                  # -gc
    maxgc = 65                  # -gc
    maxrpt = 6                  # -mr
    mt_primer3 = False          # -tm3 Use primer3-py to compute Tm?
    heterodimers = False        # Check for heterodimers formation (only if primer3-py is available)
    pr3_tm = 50                 # Tm for hairpin / heterodimer computation
    pr3_ds = -9                 # ds for hairpin / heterodimer computation

    maxover = 0.15              # -mo Amplicon size can increase by this much
    maxunder = 0.05             # -mu Amplicon size can decrease by this much
    sizemethod = "s"            # -S Method to weight amplicon size: "s" (size), "t" (targets)
    
    # Computed
    regsize = 500               # Desired amplicon size
    regmax = 0                  # Maximum amplicon size
    regmin = 0                  # Minimum amplicon size
    targetA = 0                 # Ideal position of A oligo
    targetB = 0                 # Ideal position of B oligo

    ncandout = 40               # Number of candidates to consider outside of target region
    ncandin = 10                # Number of candidates to consider inside of target region
    maxtriples = 100            # Maximum number of triples to keep for each target (sorted by best score)
    nrounds = 1000000           # -n Rounds of global optimization

    weightmt = 1.0              # -wm Weight of MT in score
    weightlen1 = 1.0            # -wl Weight of region length in score (if larger than regsize)
    weightlen2 = 1.0            # -ws Weight of region length in score (if smaller than regsize)

    def __init__(self):
        self.genes = []
        self.sequences = []

    def parseArgs(self, args):
        if "-h" in args or "--help" in args:
            return self.usage(args)
        if not checkBedtools():
            sys.stderr.write("Error: bedtools not found. This program requires bedtools to run.\n")
            sys.exit(2)

        prev = ""
        for a in args:
            if prev == "-o":
                self.outfile = a
                prev = ""
            elif prev == "-g":
                self.genelist = a
                prev = ""
            elif prev == "-u":
                self.upstream = int(a)
                prev = ""
            elif prev == "-d":
                self.downstream = int(a)
                prev = ""
            elif prev == "-s":
                self.field = int(a)
                prev = ""
            elif prev == "-tl":
                self.minmt = float(a)
                prev = ""
            elif prev == "-th":
                self.maxmt = float(a)
                prev = ""
            elif prev == "-gc":
                if "," in a:
                    parts = a.split(",")
                    if parts[0]:
                        self.mingc = float(parts[0])
                    if parts[1]:
                        self.maxgc = float(parts[1])
                prev = ""
            elif prev == "-ol":
                if "," in a:
                    parts = a.split(",")
                    if parts[0]:
                        self.minlength = int(parts[0])
                    if parts[1]:
                        self.maxlength = int(parts[1])
                else:
                    v = int(a)
                    self.minlength = v
                    self.maxlength = v
                prev = ""
            elif prev == "-wm":
                self.weightmt = float(a)
                prev = ""
            elif prev == "-wl":
                self.weightlen1 = float(a)
                prev = ""
            elif prev == "-ws":
                self.weightlen2 = float(a)
                prev = ""
            elif prev == "-w":
                self.weightlen1 = float(a)
                self.weightlen2 = float(a)
                prev = ""
            elif prev == "-n":
                self.nrounds = int(a)
                prev = ""
            elif prev == "-F":
                self.toFasta = a
                prev = ""
            elif prev == "-mo":
                self.maxover = int(a) / 100.0
                prev = ""
            elif prev == "-mu":
                self.maxunder = int(a) / 100.0
                prev = ""
            elif prev == "-S":
                self.sizemethod = a
                prev = ""
            elif prev == "-pt":
                self.pr3_tm = int(a)
                prev = ""
            elif prev == "-pd":
                self.pr3_ds = int(a)
                prev = ""
            elif prev == "-r":
                if a == "-":
                    self.reportfile = "/dev/null"
                else:
                    self.reportfile = a
                prev = ""
            elif prev == "-O":
                self.oligoFasta = a
                prev = ""
            elif prev == "-mr":
                self.maxrpt = int(a)
                prev = ""
            elif a in ["-o", "-g", "-u", "-d", "-s", "-wm", "-wl", "-ws", "-w", "-n", "-F", "-mo", "-mu", "-S", "-pt", "-pd", "-tl", "-th", "-r", "-O", "-gc", "-mr"]:
                prev = a
            elif a == "-D":
                self.heterodimers = True
            elif a == "-f":
                self.toFasta = True
            elif a == "-l":
                self.local = True
            elif a == "-L":
                self.transcripts = "longest"
            elif a == "-c":
                self.fixChroms = True
            elif a == "-nt":
                self.includeTSS = False
            elif a == "-tm3":
                if HAS_PRIMER3:
                    self.mt_primer3 = True
                else:
                    sys.stderr.write("Warning: the primer3-py package is not installed. Falling back to built-in Tm calculation.\n")
            elif self.reference is None:
                if os.path.isfile(a):
                    self.reference = a
                else:
                    sys.stderr.write("Error: genome reference file `{}' not found.\n".format(a))
                    return False
            else:
                if a[0] == '@':
                    with open(a[1:], "r") as f:
                        self.genes += [ g.strip() for g in f.read().split("\n") if g ]
                else:
                    self.genes.append(a)

        self.genes.sort()

        if not self.reference:
            return self.usage()

        if self.genelist:
            if os.path.isfile(self.genelist):
                path = self.genelist
                ext = os.path.splitext(self.genelist)[1]
            else:
                sys.stderr.write("Error: gene database file `{}' not found.\n".format(self.genelist))
                return False
        else:
            (prefix, ext, path) = self.findGeneDB()
            if not path:
                sys.stderr.write("Error: no gene database found with the name {} and one of the\n  extensions .gtf, .GTF, .gff, .gff3, .GFF, or .GFF3. Use -g to specify it.\n".format(prefix))
                return False

        if ext in [".gtf", ".GTF"]:
            self.genelist = GTFParser(path)
        elif ext in [".gff", ".gff3", ".GFF", ".GFF3"]:
            self.genelist = GFFParser(path)
        elif ext == ".txt":
            self.genelist = refFlatParser(path)
        elif ext == ".bed":
            self.genelist = BEDfileParser(path)
            self.bedmode = True
            # Change defaults for -u and -d to 250, unless user has specified one of them.
            if "-u" in args or "-d" in args:
                pass
            else:
                self.upstream = 250
                self.downstream = 250

        if self.genelist:
            self.genelist.fixChroms = self.fixChroms

        if self.genes or self.bedmode:
            return True
        else:
            return self.usage()

    def findGeneDB(self):
        base = os.path.splitext(self.reference)[0]
        name = base.split("/")[-1]
        for ext in [".gtf", ".GTF", ".gff", ".gff3", ".GFF", ".GFF3"]:
            path = base + ext
            if os.path.isfile(path):
                return (name, ext, path)
        return (name, None, None)

    def banner(self, out):
        out.write("""\x1b[1;36m***************************************************
* fengc.py - design primers for FenGC experiments *
***************************************************\x1b[0m
""")

    def usage(self, args=[]):
        self.banner(sys.stdout)
        if "weight" in args:
            sys.stdout.write("""\x1b[1mOptimization weights:\x1b[0m

This program tests all possible pairs of candidate oligos at the 5' and 3'
end of the target sequence, assigning each pair a score, and selects the
pair with the optimal score. The score is the sum of five components:

  A. Squared difference between Tm for Oligo 1 and 65.
  B. Squared difference between Tm for Oligo 2 and 65.
  C. Squared difference between Tm for Oligo 3 and 65.
  D. Square of the difference between the highest and lowest Tms.
  E. Amplicon size factor.

The first four components are multiplied by weight `wm'. Component E can be computed 
in two different ways, depending on the value of the -S argument. If the value is 
"s" (the default), E is the predicted amplicon size (i.e. the distance between the 
start of Oligo1 and the start of Oligo2) multiplied by `wl' if it is larger than 
the desired amplicon size, or by `ws' if it is smaller. 

If -S is "t", component E is the square of the distance between the start of Oligo1 
and the desired amplicon start, multipled by `wl', plus the square of the distance 
between the start of Oligo2 and the desired amplicon end, multiplied by `ws'. Use
"-h size" to view an example.

The values of the five components are added together to generate the final score.

Higher weights mean that the corresponding component has more impact on the choice
of the optimal pair, while a weight of 0 means that the corresponding component
has no impact. Examples:

* To select oligos based only on the size of the amplified region, ignoring the MT,
  use `-wm 0'.

* To avoid producing an amplified region shorter than the desired one, use `-ws 1000'
  (or any other very large number).

""")
        elif "size" in args:
            sys.stdout.write("""\x1b[1mAmplicon size weighting example\x1b[0m

In this example, the gene TSS is at position 1,000. The desired amplicon extends
from 400bp upstream of the TSS to 100b downstream, for a total size of 500bp.
Therefore, t1 is at position 600 and t2 at position 1100.


         t1                                   TSS       t2
         |              400bp                 |  100bp  |
---------[####################################+#########]--------
       ^                                               ^
       A                                               B

Imagine that OligoA is at position 570, and OligoB at position 1090, with wl=2
and ws=3. Component E of the weight is computed in the following ways.

If -S is "s":

  Predicted amplicon size is 520, which is larger than the desired size (500),
  therefore E = 520*2 = 1040.

If -S is "t":

  The distance between A and t1 is 30, and between B and t2 is 10. Therefore
  E = (30^2)*2 + (10^2)* 2 = 900*2 + 100*3 = 1800 + 300 = 2100.

""")
        elif "design" in args:
            sys.stdout.write("""\x1b[1mAmplicon design:\x1b[0m

These are the steps used by FenGC to design an amplicon starting from the
TSS position of a gene.

* The `desired amplicon' by default extends from 400bp upstream of the TSS
  to 100bp downstream of the TSS. These values can be changed with the -u
  and -d options respectively.

* The two options -mo and -mu represents how much the amplicon size is allowed
  to grow or shrink, in percentage. The default values are 15 and 5, which means
  that the amplicon size can range from 475bp to 575bp.

* The program looks for all potential oligos around the left and right boundaries
  of the desired amplicon, and ranks them according to the procedure described in
  `-h weight'. This is done separately for OligoA, OligoB, and OligoC, and the
  set with the best score is chosen as the optimal triple (this can be changed
  later if global optimization is enabled).

* The program forces OligoA to be upstream of the TSS and OligoB to be downstream
  of it. To remove this constraing (and allow amplicons that do not cover the 
  TSS) add the -nt option.

""")              

        else:
            sys.stdout.write("""
\x1b[1mUsage: fengc.py [options] reference genes...\x1b[0m

where `reference' is a FASTA file containing the genome reference, and `genes' is one or
more gene names, or (if preceded by @) a file containing gene names, one per line.

\x1b[1mGeneral options:\x1b[0m

  -g G  | Use gene database G (in GTF or GFF format). Default: a file with the same
          name as the reference file, with extension .gtf or .gff (lower or uppercase).
  -c    | Add "chr" at the beginning of chromosome names that don't have it.
  -L    | For genes having multiple isoforms, select longest one (default: select all).
  -o O  | Write output to file O (default: standard output).
  -f    | Write target sequences to separate FASTA files.
  -F F  | Write target sequences to single FASTA file F.
  -O O  | Write oligo sequences to a FASTA file named O.
  -r R  | Write report to file R (default: {}).
  -D    | Check for potential heterodimers.

\x1b[1mTm options:\x1b[0m
  
  -tl L   | Set minimum allowed Tm (default: {}).
  -th H   | Set maximum allowed Tm (default: {}).
  -ol L,H | Choose oligos in the size range L - H (default {},{}).
  -gc L,H | Set GC range to L - H (default: {},{}).
  -tm3    | Use primer3 method to compute Tm (default: builtin method).
  -pt T   | Limit temperature for hairpin/heterodimers (default: {}).
  -pd D   | Limit Ds for hairpin/heterodimers (default: {}).
  -mr R   | Maximum number of repeat bases (default: {})."

\x1b[1mDesign options (see -h design and -h size for details):\x1b[0m

  -u U  | Number of bp upstream of TSS (default: {}).
  -d D  | Number of bp downstream of TSS (default: {}).
  -s S  | Number of bp upstream/downstream of regions of interest (default: {}).
  -mo O | Maximum % increase of amplicon size (default: {}%).
  -mu U | Maximum % decrease of amplicon size (default: {}%).
  -S A  | Amplicon sizing method, one of `s' or `t' (default: {}).
  -nt   | Do NOT force amplicon to contain TSS.

\x1b[1mWeight options (see -h weight for details):\x1b[0m

  -wm W | Set weight for Tm penalty (default: {}).
  -wl W | Set weight for region length when larger than target (default: {}).
  -ws W | Set weight for region length when smaller than target (default: {}).
  -w  W | Set both -wl and -ws to W.

\x1b[1mOptimization options:\x1b[0m

  -n N | Perform N rounds of global optimization (default: {}).
  -l   | Local only - do not perform global optimization.

\x1b[1mMore details:\x1b[0m

  fengc -h design
  fengc -h size
  fengc -h weight

\x1b[1mRequirements:\x1b[0m

  bedtools
  primer3-py (only if using the -tm3 or -D options).
    Available at: https://pypi.org/project/primer3-py/

(c) 2020, University of Florida Research Foundation.

""".format(self.reportfile, self.minmt, self.maxmt, self.minlength, self.maxlength, self.mingc, self.maxgc, self.pr3_tm, self.pr3_ds, self.maxrpt, self.upstream, self.downstream, self.field, 
           int(self.maxover * 100), int(self.maxunder * 100), self.sizemethod,
           self.weightmt, self.weightlen1, self.weightlen2, 
           self.nrounds))
        return False

    def run(self):
        with open(self.reportfile, "w") as rout:
            rep = Report(rout)
            self.gnames = []
            self.badgenes = []
            self.gcoords = []
            self.regsize  = self.upstream + self.downstream # Size of target region
            self.regmax = self.regsize + int(self.regsize * self.maxover)
            self.regmin = self.regsize - int(self.regsize * self.maxunder)
            self.targetA = self.field
            self.targetB = self.field + self.regsize

            if self.bedmode:
                sys.stderr.write("BED mode                    True\n")
            else:
                sys.stderr.write("Input genes:                {}\n".format(len(self.genes)))
                sys.stderr.write("Transcript selection:       {}\n".format(self.transcripts))
            sys.stderr.write("Amplicon size range:        {} - {}\n".format(self.regmin, self.regmax))
            if self.bedmode:
                sys.stderr.write("Optimal amplicon positions: CTR-{}, CTR+{}\n".format(self.upstream, self.downstream))
            else:
                sys.stderr.write("Optimal amplicon positions: TSS-{}, TSS+{}\n".format(self.upstream, self.downstream))
            sys.stderr.write("Tm range:                   {} - {}\n".format(self.minmt, self.maxmt))
            sys.stderr.write("GC% range:                  {} - {}\n".format(self.mingc, self.maxgc))
            sys.stderr.write("Global optimization:        {}\n".format("disabled" if self.local else "enabled"))
            if not self.local:
                sys.stderr.write("  rounds:                   {}\n".format(self.nrounds))
            sys.stderr.write("Heterodimers check:         {}\n".format("enabled" if self.heterodimers else "disabled"))
            sys.stderr.write("Output files:\n")
            sys.stderr.write("  Oligo design table:       {}\n".format(self.outfile))
            sys.stderr.write("  Report file:              {}\n".format(self.reportfile))
            if self.toFasta:
                sys.stderr.write("  Amplified sequences:      {}\n".format(self.toFasta))
            if self.oligoFasta:
                sys.stderr.write("  Oligo sequences:          {}\n".format(self.oligoFasta))
            sys.stderr.write("\n")
            self.seqman = SequenceManager(self.reference)
            self.genelist.load()

            sys.stderr.write("\n\n\x1b[1mLoading target sequences:\x1b[0m\n")
            self.findGeneCoords()
            if self.badgenes:
                rep.badGenes(self.badgenes)
            #sys.stderr.write("{}\n".format(self.gcoords))
            self.sequences = self.seqman.getSequences(self.gcoords) # This will be a list of Sequence objects

            sys.stderr.write("\n\x1b[1mFinding oligo triples - best triple for each sequence:\x1b[0m\n")
            sys.stderr.write("  \x1b[4mGene\x1b[0m                        \x1b[4mSize\x1b[0m    \x1b[4m Oligo 1 \x1b[0m   \x1b[4m Oligo 2 \x1b[0m   \x1b[4m Oligo 3 \x1b[0m\n")
            sys.stderr.write("                                      MT    %GC   MT    %GC   MT    %GC\n")
            regstart = self.field                # Start of target region
            regend   = self.field + self.regsize # End of target region

            badoligos = []
            for i in range(len(self.sequences)):
                seq = self.sequences[i]
                if not seq.seq: # Genes for which we don't have a reference sequence are skipped immediately
                    seq.valid = False
                    sys.stderr.write("  {:20}  -- no valid oligos found --\n".format(seq.name))
                    badoligos.append(seq.name)
                    continue

                seq.name = self.gnames[i]
                if seq.strand == '-':
                    seq.seq = revcomp(seq.seq)

                # Here's where we find the best oligo triple for this sequence
                seq.triples = self.findOptimalOligos(seq, regstart, regend) # List of triples, best first

                if seq.triples:
                    best = seq.best = seq.triples[0]
                    sys.stderr.write("  {:28}{}bp   {:.1f}  {}%   {:.1f}  {}%   {:.1f}  {}%\n".format(
                        seq.name, best.size,
                        best.oligo1.mt, int(100*best.oligo1.gcperc), 
                        best.oligo2.mt, int(100*best.oligo2.gcperc), 
                        best.oligo3.mt, int(100*best.oligo3.gcperc)))
                    #offset = seq.start
                    #sys.stderr.write("{} - {} - {}\n".format(best.oligo1.start+offset, best.oligo2.start+offset, best.oligo3.start+offset))
                else:
                    seq.valid = False
                    sys.stderr.write("  {:20}  -- no valid oligos found --\n".format(seq.name))
                    badoligos.append(seq.name)

            if badoligos:
                rep.badOligos(badoligos)

            # Global optimization
            if not self.local:
                self.globalOptimization()

            # Check for heterodimers
            if HAS_PRIMER3 and self.heterodimers:
                self.checkHeterodimers(rep)

            # Writing output
            sys.stderr.write("\n\x1b[1mWriting output:\x1b[0m\n")
            nbad = self.writeOutput()
            if self.toFasta:
                self.writeFastas()
            if self.oligoFasta:
                self.writeOligoFastas()
            sys.stderr.write("\n\x1b[1mSummary:\x1b[0m\n")
            sys.stderr.write("  {} input transcripts\n".format(len(self.sequences)))
            sys.stderr.write("  Oligo design \x1b[32msuccessful\x1b[0m for \x1b[1m{}\x1b[0m transcripts\n".format(len(self.sequences) - nbad))
            sys.stderr.write("  Oligo design \x1b[31mfailed    \x1b[0m for \x1b[1m{}\x1b[0m transcripts\n".format(nbad))
            sys.stderr.write("\n")

    def findGeneCoords(self):
        self.badgenes = []

        if self.bedmode:
            for gname in self.genelist.genes.keys():
                tx = self.genelist.get(gname).transcripts[0]
                midpoint = int((tx[2] + tx[3]) / 2)
                self.gcoords.append([tx[1], midpoint - self.upstream-self.field, midpoint + self.downstream + self.field, tx[4]])
                self.gnames.append(tx[0])
                sys.stderr.write("  {:30} {}:{}-{}:{}\n".format(tx[0], tx[1], midpoint-self.upstream, midpoint+self.downstream, tx[4]))
        for g in self.genes:
            if ":" in g:
                gname = g.split(":")[0]
                acc = g
            else:
                gname = g
                acc = None
            gene = self.genelist.get(gname)
            if gene:
                maxlength = 0
                wanted = []
                for tx in gene.transcripts:
                    if acc:                           # If we're interested in a single transcript
                        if tx[0] == acc:              # and it's this one
                            wanted.append(tx)         # we're done.
                            break
                    elif self.transcripts == "all":   # If we want all of them
                        wanted.append(tx)             # add this one
                    else:
                        length = tx[3] - tx[2]        # If we want the longest one
                        if length > maxlength:        # and this one is longest
                            wanted = [tx]             # replace previous one with this
                            maxlength = tx[3] - tx[2] # and set new maxlength

                for tx in wanted:
                    #sys.stderr.write("{}: TSS = {}\n".format(tx[0], tx[2] if tx[4] == "+" else tx[3]))
                    self.gnames.append(tx[0])
                    if tx[4] == '+':
                        self.gcoords.append([tx[1], tx[2] - self.upstream - self.field, tx[3] + self.downstream + self.field, tx[4]])
                    else:
                        self.gcoords.append([tx[1], tx[2] - self.downstream - self.field, tx[3] + self.upstream + self.field, tx[4]])
                    sys.stderr.write("  {:30} {}:{}-{}:{}\n".format(tx[0], tx[1], tx[2], tx[3], tx[4]))
            else:
                self.badgenes.append(g)

    def checkHeterodimers(self, rep):
        sys.stderr.write("\n\x1b[1mChecking for potential heterodimers:\x1b[0m\n")
        hets = []
        l = len(self.sequences)
        ntests = l * (l - 1) / 2
        ndone = 0
        for i in range(l):
            g1  = self.sequences[i].name
            one = self.sequences[i].best
            if one:
                self.checkHeterodimersAux(one.oligo1, one.oligo2, g1, g1, hets)
                self.checkHeterodimersAux(one.oligo1, one.oligo3, g1, g1, hets)
                self.checkHeterodimersAux(one.oligo2, one.oligo3, g1, g1, hets)
                for j in range(i+1, l):
                    g2  = self.sequences[j].name
                    two = self.sequences[j].best
                    if two:
                        self.checkHeterodimersAux(one.oligo1, two.oligo1, g1, g2, hets)
                        self.checkHeterodimersAux(one.oligo1, two.oligo2, g1, g2, hets)
                        self.checkHeterodimersAux(one.oligo1, two.oligo3, g1, g2, hets)
                        self.checkHeterodimersAux(one.oligo2, two.oligo1, g1, g2, hets)
                        self.checkHeterodimersAux(one.oligo2, two.oligo2, g1, g2, hets)
                        self.checkHeterodimersAux(one.oligo2, two.oligo3, g1, g2, hets)
                        self.checkHeterodimersAux(one.oligo3, two.oligo1, g1, g2, hets)
                        self.checkHeterodimersAux(one.oligo3, two.oligo2, g1, g2, hets)
                        self.checkHeterodimersAux(one.oligo3, two.oligo3, g1, g2, hets)
                    ndone += 1
                    sys.stderr.write("\x1b[G  Progress: {}%".format(int(100.0 * ndone / ntests)))
        sys.stderr.write("\x1b[G  Done - {} potential heterodimers found.\n".format(len(hets)))
        if hets:
            rep.heterodimers(hets)

    def checkHeterodimersAux(self, o1, o2, g1, g2, hets):
        het_info = primer3.calcHeterodimer(o1.sequence, o2.sequence)
        if het_info.structure_found == True:
            #sys.stderr.write("Heterodimer: {} {}\n".format(het_info.ds, het_info.tm))
            if het_info.ds <= self.pr3_ds and het_info.tm >= self.pr3_tm:
                hets.append( (g1, o1.sequence, g2, o2.sequence, het_info.tm) )
                if "H" not in o1.flags:
                    o1.flags += "H"
                if "H" not in o2.flags:
                    o1.flags += "H"

    def writeOutput(self):
        nbad = 0
        sys.stderr.write("  Writing oligo information for {} genes to {}.\n".format(len(self.sequences), self.outfile))
        with open(self.outfile, "w") as out:
            out.write("Gene\tOrientation\tGenomic coordinates\tOligo pos (from TSS)\tPCR Product Size\tPCR Product GC%\tO1-sequence\tO1-length\tO1-MT\tO1-GC%\tO2-sequence\tO2-length\tO2-MT\tO2-GC%\tO3-sequence\tO3-length\tO3-MT\tO3-GC%\n")
            for i in range(len(self.sequences)):
                seq = self.sequences[i]
                if not seq.valid:
                    nbad += 1
                    continue
                best = seq.best
                #sys.stderr.write("{} - {} - {}\n".format(best.oligo1.start, best.oligo2.start, best.oligo3.start))
                out.write("{}\t{}\t{}:{:,}-{:,}\t{} {}\t{}\t{}\t{}\t{}\t{:.1f}\t{}\t{}\t{}\t{:.1f}\t{}\t{}\t{}\t{:.1f}\t{}\n".format(
                    seq.name, "F" if seq.strand == "+" else "R", seq.chrom, seq.start, seq.end, 
                    best.oligo1.start+1-self.field-self.upstream, best.oligo2.start+1-self.field-self.upstream,
                    best.size, int(100*seq.bestGC()),
                    revcomp(best.oligo1.sequence)+U1, len(best.oligo1.sequence), best.oligo1.mt, int(100*best.oligo1.gcperc),
                    revcomp(best.oligo2.sequence)+U1, len(best.oligo2.sequence), best.oligo2.mt, int(100*best.oligo2.gcperc),
                    U2+revcomp(best.oligo3.sequence), len(best.oligo3.sequence), best.oligo3.mt, int(100*best.oligo3.gcperc)))
        return nbad

    def writeFastas(self):
        if self.toFasta is True:
            for seq in self.sequences:
                if seq.best:
                    seq.write(seq.name + ".fa", label=seq.name + " " + seq.coordinates(), start=seq.best.oligo1.start+1, end=seq.best.oligo3.end)
                    #seq.write(seq.name + ".fa", label=seq.name + " " + seq.coordinates(), start=seq.best.oligo1.start, end=seq.best.oligo2.end)
        else:
            with open(self.toFasta, "w") as out:
                for seq in self.sequences:
                    if seq.best:
                        #sys.stderr.write("{}-{}\n".format(seq.best.oligo1.start+1, seq.best.oligo3.end))
                        seq.writes(out, label=seq.name + " " + seq.coordinates(), start=seq.best.oligo1.start+1, end=seq.best.oligo3.end)
                        #seq.writes(out, label=seq.name + " " + seq.coordinates(), start=seq.best.oligo1.start, end=seq.best.oligo2.end)

    def writeOligoFastas(self):
        with open(self.oligoFasta, "w") as out:
            for seq in self.sequences:
                best = seq.best
                if best:
                    out.write(">{}_A\n{}\n".format(seq.name, best.oligo1.sequence))
                    out.write(">{}_B\n{}\n".format(seq.name, best.oligo2.sequence))
                    out.write(">{}_C\n{}\n".format(seq.name, best.oligo3.sequence))

    # Find optimal oligo length in a sequence giving MT as close as possible to 65 degrees
    # Returns an Oligo object

    def findOptimalMT(self, seq, start, direction):
        # print (start, direction, minlen, maxlen)
        bestlen = 0
        bestmt = 100
        delta = 100
        for l in range(self.minlength, self.maxlength):
            if direction == 1:
                oligo = seq[start:start+l] ## TODO: optimize - we just need to look at the base being added
            else:
                oligo = seq[start-l+1:start+1]
            gc = ngc(oligo)
            gcperc = 100.0 * gc / l

            # We want GC% between 40 and 65
            if gcperc < self.mingc or gcperc > self.maxgc:
                continue

            # Count lower-case (i.e. repeatmasked) bases
            nl = nlower(oligo)
            if nl >= self.maxrpt:
                continue

            # Check hairpin formation
            if self.mt_primer3:
                hairpin_info = primer3.calcHairpin(oligo)
                if hairpin_info.structure_found == True:
                    if hairpin_info.ds <= self.pr3_ds and hairpin_info.tm >= self.pr3_tm:
                        continue

            # Compute MT
            if self.mt_primer3:
                thismt = primer3.calcTm(oligo)
            else:
                thismt = mt(gc, l)

            # If we're closer to 65 than the current best, keep this
            if abs(thismt-65) < delta:
                bestlen = l
                bestmt = thismt
                delta = abs(thismt-65)

        if self.minmt <= bestmt <= self.maxmt:
            return Oligo(seq, start, direction, bestlen, bestmt)
        else:
            return None

    def findOptimalOligos(self, seq, regstart, regend):
        tsspos = regstart + self.upstream
        triples = []

        positions = seq.findVT()
        #sys.stderr.write("{}\n".format((regstart, regend, tsspos)))
        pos1 = findCandidatePositions(positions, self.ncandout, regstart, -1) + findCandidatePositions(positions, self.ncandin, regstart, 1)
        pos1.sort(key=lambda p: abs(regstart - p))
        #sys.stderr.write("{}\n".format(pos1))
        pos2 = findCandidatePositions(positions, self.ncandout, regend, 1) + findCandidatePositions(positions, self.ncandin, regend, -1)
        pos2.sort(key=lambda p: abs(regend - p))
        if self.includeTSS:
            pos1 = [ p for p in pos1 if p < tsspos ]
            pos2 = [ p for p in pos2 if p > tsspos ]


        #sys.stderr.write("{}\n".format(pos2))
        oligos1 = [ self.findOptimalMT(seq.seq, start, 1) for start in pos1 ]
        oligos2 = [ self.findOptimalMT(seq.seq, start, 1) for start in pos2 ]

        # Find best pair
        best = None
        bestscore = 1e6
        for a in oligos1:
            for b in oligos2:
                if a and b:
                    # See if we can find Oligo3
                    c = self.findOptimalMT(seq.seq, b.start, -1)
                    if c:
                        t = Triple(a, b, c)
                        if self.regmin <= t.size <= self.regmax:
                            t.score = self.score(t)
                            triples.append(t)
        triples.sort(key=lambda t: t.score)
        return triples[:self.maxtriples]

    def score(self, t):
        s1 = (t.oligo1.mt - 65)**2
        s2 = (t.oligo2.mt - 65)**2
        s3 = (t.oligo3.mt - 65)**2
        s4 = t.tmspread**2
        if self.sizemethod == "s":
            size = t.size - self.regsize
            if size > 0:
                s5 = self.weightlen1 * size
            else:
                s5 = self.weightlen2 * size
        elif self.sizemethod == "t":
            s5 = self.weightlen1 * (t.oligo1.start - self.targetA)**2 + self.weightlen2 * (t.oligo2.start - self.targetB)**2
        return self.weightmt*(s1 + s2 + s3 + s4) + s5

    def randomTripleset(self):
        """Pick a triple at random for each sequence (storing it into rndtriple)."""
        for seq in self.sequences:
            if seq.valid:
                seq.randomTriple()

    def triplesetVariance(self):
        """Return the variance of the Tms of the current set of random triples."""
        tms = []
        for seq in self.sequences:
            if seq.valid:
                t = seq.rndtriple
                tms.append(t.oligo1.mt)
                tms.append(t.oligo2.mt)
                tms.append(t.oligo3.mt)
        return var(tms)

    def setBestTripleset(self):
        for seq in self.sequences:
            if seq.valid:
                seq.best = seq.rndtriple

    def globalOptimization(self):
        sys.stderr.write("\n\x1b[1mGlobal optimization:\x1b[0m\n")
        sys.stderr.write("  Performing MonteCarlo optimization, {:,} rounds (ctrl-c to interrupt)\n".format(self.nrounds))
        sys.stderr.flush()
        bestvar = 100000
        try:
            for i in range(self.nrounds):
                update = False
                self.randomTripleset()
                var = self.triplesetVariance()
                if var < bestvar:
                    self.setBestTripleset()
                    bestvar = var
                    update = True
                if update or i % 10000 == 0:
                    sys.stderr.write("\x1b[G  Round: {:,}, Best variance: {:.3f}".format(i, bestvar))
                    sys.stderr.flush()
        except KeyboardInterrupt:
            pass
        sys.stderr.write("\n  {:,} rounds done, best variance={:.3f}\n".format(i+1, bestvar))

### Let's get things started

if __name__ == "__main__":
    args = sys.argv[1:]
    M = Main()
    if M.parseArgs(args):
        M.banner(sys.stderr)
        M.run()

### Test with:

## ./fengc.py -f /ufrc/data/reference/icbr/GRCh38/GRCh38.fa -g /ufrc/data/reference/icbr/GRCh38/Homo_sapiens.GRCh38.95.gtf APOE KDM6A TLR4
