#!/usr/bin/env python

import sys
import csv
import os.path
import random
from numpy import var
import subprocess as sp
import importlib

try:
    importlib.import_module('primer3')
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

BASES = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}

def revcomp(s):
    bases = [ BASES[b] for b in s[::-1] ]
    return "".join(bases)

# Number of GC bases in a sequence

def ngc(seq):
    return seq.count("C") + seq.count("G")

# Find candidate positions to test

def findCandidatePositions(positions, howmany, start, direction):
    """Return `howmany' positions from the list of positions starting at `start' in direction `direction'
(1 = left to right, -1 = right to left)."""
    
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

    def write(self, filename, label=None):
        with open(filename, "w") as out:
            self.writes(out, label=label)

    def writes(self, out, label=None):
        out.write(">{}\n{}\n".format(label or self.name, self.seq))

    def findVT(self):
        positions = []
        for i in range(1, len(self.seq)):
            if self.seq[i] == 'T' and self.seq[i-1] != 'T':
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
        seqs = []
        proc = sp.Popen("bedtools getfasta -fi {} -bed /dev/stdin".format(self.reference), shell=True, stdin=sp.PIPE, stdout=sp.PIPE)
        for intv in intervals:
            proc.stdin.write("{}\t{}\t{}\n".format(intv[0], intv[1], intv[2]))
        proc.stdin.close()
        for intv in intervals:
            name = proc.stdout.readline().rstrip("\n")
            seq = proc.stdout.readline().rstrip("\n").upper()
            s = Sequence(name, seq)
            s.chrom = intv[0]
            s.start = intv[1]
            s.end   = intv[2]
            s.strand = intv[3]
            seqs.append(s)
        return seqs

# Simple GTF parser

class GTFParser(object):
    filename = ""
    genes = {}

    def __init__(self, filename):
        self.filename = filename
        self.genes = {}

    def fixChrom(self, c):
        if c.startswith("chr"):
            return c
        else:
            return "chr" + c

    def getName(self, annots):
        parts = annots.split(";")
        for p in parts:
            p = p.lstrip()
            if p.startswith("gene_name"):
                return p[10:].strip('"')
        return None

    def load(self):
        sys.stderr.write("Loading genes from GTF file {}: \x1b[s".format(self.filename))
        n = 0
        with open(self.filename, "r") as f:
            c = csv.reader(f, delimiter='\t')
            for line in c:
                if line[0][0] == '#':
                    continue
                if line[2] != "gene":
                    continue
                chrom = self.fixChrom(line[0])
                start = int(line[3])
                end = int(line[4])
                strand = line[6]
                name = self.getName(line[8])
                if name:
                    self.genes[name] = [chrom, start, end, strand]
                    n += 1
                    if n % 5000 == 0:
                        sys.stderr.write("\x1b[u{}".format(n))
        sys.stderr.write("\x1b[u{} genes loaded.\n".format(len(self.genes)))

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

## Main

class Main(object):
    reference = None            # Name of reference file
    genelist = None             # -g Name of genes DB
    genes = []
    sequences = []              # List of Sequence objects
    outfile = "/dev/stdout"     # -o
    seqman = None
    toFasta = False             # -F, -f
    local = False               # -l If true, skips global optimization.
    
    upstream = 400              # -u
    downstream = 100            # -d
    field = 2000                # -s
    minlength = 8
    maxlength = 30
    minmt = 62
    maxmt = 68
    mt_primer3 = False          # -tm3 Use primer3-py to compute Tm?

    maxover = 0.15              # -mo Amplicon size can increase by this much
    maxunder = 0.05             # -mu Amplicon size can decrease by this much
    sizemethod = "s"            # -S Method to weight amplicon size: "s" (size), "t" (targets)

    # Computed
    regsize = 500               # Desired amplicong size
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
                self.maxover = 100.0 / int(a)
                prev = ""
            elif prev == "-mu":
                self.maxunder = 1,0 - (100.0 / int(a))
                prev = ""
            elif prev == "-S":
                self.sizemethod = a
                prev = ""
            elif a in ["-o", "-g", "-u", "-d", "-s", "-wm", "-wl", "-ws", "-w", "-n", "-F", "-mo", "-mu", "-S"]:
                prev = a
            elif a == "-f":
                self.toFasta = True
            elif a == "-l":
                self.local = True
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
                        self.genes += [ g.strip() for g in f.read().split("\n") ]
                else:
                    self.genes.append(a)

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

        if self.genes:
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

The first three components are multiplied by weight `wm'. Component E can be computed 
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
         [              400bp                 |  100bp  ]
----------------------------------------------+------------------
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
        else:
            sys.stdout.write("""
\x1b[1mUsage: fengc.py [options] reference genes...\x1b[0m

where `reference' is a FASTA file containing the genome reference, and `genes' is one or
more gene names, or (if preceded by @) a file containing gene names, one per line.

\x1b[1mGeneral options:\x1b[0m

  -g G | Use gene database G (in GTF or GFF format). Default: a file with the same
         name as the reference file, with extension .gtf or .gff (lower or uppercase).
  -o O | Write output to file O (default: standard output).
  -f   | Write target sequences to separate FASTA files.
  -F F | Write target sequences to single FASTA file F.
  -tm3 | Use primer3 method to compute Tm (default: builtin method).

\x1b[1mDesign options (see -h design for details):\x1b[0m

  -u U  | Number of bp upstream of TSS (default: {}).
  -d D  | Number of bp downstream of TSS (default: {}).
  -s S  | Number of bp upstream/downstream of regions of interest (default: {}).
  -mo O | Maximum % increase of amplicon size (default: {}%).
  -mu U | Maximum % decrease of amplicon size (default: {}%).
  -S A  | Amplicon sizing method, one of `s' or `t' (default: {}).

\x1b[1mWeight options (see -h weight for details):\x1b[0m

  -wm W | Set weight for Tm penalty (default: {}).
  -wl W | Set weight for region length when larger than target (default: {}).
  -ws W | Set weight for region length when smaller than target (default: {}).
  -w  W | Set both -wl and -ws to W.

\x1b[1mOptimization options:\x1b[0m

  -n N | Perform N rounds of global optimization (default: {}).
  -l   | Local only - do not perform global optimization.

\x1b[1mRequirements:\x1b[0m

  bedtools

  primer3-py (only if useing the -tm3 option).
    Available at: https://pypi.org/project/primer3-py/

(c) 2020, University of Florida Research Foundation.

""".format(self.upstream, self.downstream, self.field, 
           int(self.maxover * 100), int(self.maxunder * 100), self.sizemethod,
           self.weightmt, self.weightlen1, self.weightlen2, 
           self.nrounds))
        return False

    def run(self):
        gnames = []
        gcoords = []
        self.regsize  = self.upstream + self.downstream # Size of target region
        self.regmax = self.regsize + int(self.regsize * self.maxover)
        self.regmin = self.regsize - int(self.regsize * self.maxunder)
        self.targetA = self.field
        self.targetB = self.field + self.regsize

        self.seqman = SequenceManager(self.reference)
        self.genelist.load()

        sys.stderr.write("\n\x1b[1mLoading target sequences:\x1b[0m\n")
        for g in self.genes:
            coords = self.genelist.get(g)
            if coords:
                gnames.append(g)
                if coords[3] == '+':
                    gcoords.append([coords[0], coords[1] - self.upstream - self.field, coords[1] + self.downstream + self.field, coords[3]])
                else:
                    gcoords.append([coords[0], coords[2] - self.downstream - self.field, coords[2] + self.upstream + self.field, coords[3]])
                sys.stderr.write("  {:20} {}:{}-{}:{}\n".format(g, coords[0], coords[1], coords[2], coords[3]))

        self.sequences = self.seqman.getSequences(gcoords) # This will be a list of Sequence objects

        sys.stderr.write("\n\x1b[1mFinding oligo triples - best triple for each sequence:\x1b[0m\n")
        sys.stderr.write("  \x1b[4mGene\x1b[0m                \x1b[4mSize\x1b[0m    \x1b[4m Oligo 1 \x1b[0m   \x1b[4m Oligo 2 \x1b[0m   \x1b[4m Oligo 3 \x1b[0m\n")
        sys.stderr.write("                              MT    %GC   MT    %GC   MT    %GC\n")
        regstart = self.field                # Start of target region
        regend   = self.field + self.regsize # End of target region
        for i in range(len(gnames)):
            seq = self.sequences[i]
            if gcoords[i][3] == '-':
                seq.seq = revcomp(seq.seq)
            seq.name = gnames[i]

            # Here's where we find the best oligo triple for this sequence
            seq.triples = self.findOptimalOligos(seq, regstart, regend) # List of triples, best first

            if seq.triples:
                best = seq.best = seq.triples[0]
                sys.stderr.write("  {:20}{}bp   {:.1f}  {}%   {:.1f}  {}%   {:.1f}  {}%\n".format(
                    gnames[i], best.size,
                    best.oligo1.mt, int(100*best.oligo1.gcperc), 
                    best.oligo2.mt, int(100*best.oligo2.gcperc), 
                    best.oligo3.mt, int(100*best.oligo3.gcperc)))
            else:
                seq.valid = False
                sys.stderr.write("  {:20}  -- no valid oligos found --\n".format(gnames[i]))

        # Global optimization
        if not self.local:
            self.globalOptimization()

        # Writing output
        sys.stderr.write("\n\x1b[1mWriting output:\x1b[0m\n")
        self.writeOutput()
        if self.toFasta:
            self.writeFastas(gnames)

    def writeOutput(self):
        sys.stderr.write("  Writing oligo information for {} genes to {}.\n".format(len(self.sequences), self.outfile))
        with open(self.outfile, "w") as out:
            out.write("Gene\tOrientation\tGenomic coordinates\tOligo positions\tPCR Product Size\tPCR Product GC%\tO1-sequence\tO1-length\tO1-MT\tO1-GC%\tO2-sequence\tO2-length\tO2-MT\tO2-GC%\tO3-sequence\tO3-length\tO3-MT\tO3-GC%\n")
            for i in range(len(self.sequences)):
                seq = self.sequences[i]
                if not seq.valid:
                    continue
                best = seq.best
                out.write("{}\t{}\t{}:{:,}-{:,}\t{}-{}\t{}\t{}\t{}\t{}\t{:.1f}\t{}\t{}\t{}\t{:.1f}\t{}\t{}\t{}\t{:.1f}\t{}\n".format(
                    seq.name, "F" if seq.strand == "+" else "R", seq.chrom, seq.start, seq.end, best.oligo1.start+1, best.oligo2.start+1,
                    best.size, int(100*seq.bestGC()),
                    revcomp(best.oligo1.sequence)+"/"+U1, len(best.oligo1.sequence), best.oligo1.mt, int(100*best.oligo1.gcperc),
                    revcomp(best.oligo2.sequence)+"/"+U1, len(best.oligo2.sequence), best.oligo2.mt, int(100*best.oligo2.gcperc),
                    revcomp(best.oligo3.sequence)+"/"+U2, len(best.oligo3.sequence), best.oligo3.mt, int(100*best.oligo3.gcperc)))

    def writeFastas(self, gnames):
        if self.toFasta is True:
            for i in range(len(gnames)):
                seq = self.sequences[i]
                seq.write(gnames[i] + ".fa", label=gnames[i] + " " + seq.coordinates())
        else:
            with open(self.toFasta, "w") as out:
                for i in range(len(gnames)):
                    seq = self.sequences[i]
                    seq.writes(out, label=gnames[i] + " " + seq.coordinates())

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
            gcperc = 1.0 * gc / l

            # We want GC% between 40 and 65
            if gcperc < 0.4 or gcperc > 0.65:
                continue

            # Check hairpin formation
            if self.mt_primer3:
                hairpin_info = primer3.calcHairpin(oligo)
                if hairpin_info.structure_found = True and hairpin_info.ds <= -9 and hairpin_info >= 50:
                    continue

            # Compute MT
            if self.mt_primer3:
                thismt = primer3.caclTm(oligo)
            else:
                thismt = mt(gc, l)

            # If we're closer to 65 than the current best, keep this
            if abs(thismt-65) < delta:
                bestlen = l
                bestmt = thismt
                delta = abs(thismt-65)

        if 60 <= bestmt < 70:
            return Oligo(seq, start, direction, bestlen, bestmt)
        else:
            return None

    def findOptimalOligos(self, seq, regstart, regend):
        triples = []
        positions = seq.findVT()
        pos1 = findCandidatePositions(positions, self.ncandout, regstart, -1) + findCandidatePositions(positions, self.ncandin, regstart, 1)
        pos1.sort(key=lambda p: abs(regstart - p))
        pos2 = findCandidatePositions(positions, self.ncandout, regend, 1) + findCandidatePositions(positions, self.ncandin, regend, -1)
        pos2.sort(key=lambda p: abs(regend - p))
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
        sys.stderr.write("  \x1b[s")
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
                    sys.stderr.write("\x1b[uRound: {:,}, Best variance: {:.3f}".format(i, bestvar))
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
