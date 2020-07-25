#!/usr/bin/env python

import sys
import csv
import os.path
import subprocess as sp

## Utilities

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
    bases = [ BASES[b.upper()] for b in s[::-1] ]
    return "".join(bases)

# Number of GC bases in a sequence

def ngc(seq):
    return seq.count("C") + seq.count("G")

# Find optimal oligo length in a sequence giving MT as close as possible to 65 degrees
# Returns an Oligo object

def findOptimalMT(seq, start, direction, minlen, maxlen):
    # print (start, direction, minlen, maxlen)
    bestlen = 0
    bestmt = 100
    delta = 100
    for l in range(minlen, maxlen):
        if direction == 1:
            oligo = seq[start:start+l] ## TODO: optimize - we just need to look at the base being added
        else:
            oligo = seq[start-l+1:start+1]
        gc = ngc(oligo)
        gcperc = 1.0 * gc / l

        # We want GC% between 40 and 65
        if gcperc < 0.4 or gcperc > 0.65:
            continue

        # Compute MT
        thismt = mt(gc, l)

        # If we're closer to 65 than the current best, keep this
        if abs(thismt-65) < delta:
            bestlen = l
            bestmt = thismt
            delta = abs(thismt-65)

    if bestmt < 100:
        return Oligo(seq, start, direction, bestlen, bestmt)
    else:
        return None

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

# Extracting sequences for regions

class Sequence(object):
    name = ""
    seq = ""
    extra = None

    def __init__(self, name, seq):
        self.name = name
        self.seq = seq

    def write(self, filename):
        with open(filename, "w") as out:
            out.write(">{}\n{}\n".format(self.name, self.seq))

    def findVT(self):
        positions = []
        for i in range(1, len(self.seq)):
            if self.seq[i] == 'T' and self.seq[i-1] != 'T':
                positions.append(i)
        return positions

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
        for _ in intervals:
            name = proc.stdout.readline().rstrip("\n")
            seq = proc.stdout.readline().rstrip("\n")
            seqs.append(Sequence(name, seq))
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
    genelist = None
    genes = []
    sequences = []              # List of Sequence objects
    outfile = "/dev/stdout"
    reference = ""
    seqman = None
    toFasta = False

    upstream = 400
    downstream = 100
    regsize = 500
    field = 2000
    minlength = 8
    maxlength = 30
    minmt = 62
    maxmt = 68

    ncandout = 40               # Number of candidates to consider outside of target region
    ncandin = 10                # Number of candidates to consider inside of target region

    weightmt = 1.0              # Weight of MT in score
    weightlen1 = 1.0            # Weight of region length in score (if larger than regsize)
    weightlen2 = 1.0            # Weight of region length in score (if smaller than regsize)

    def __init__(self):
        self.genes = []
        self.sequences = []

    def parseArgs(self, args):
        if "-h" in args or "--help" in args:
            return False
        prev = ""
        for a in args:
            if prev == "-o":
                self.outfile = a
                prev = ""
            elif prev == "-r":
                self.reference = a
                prev = ""
            elif prev == "-u":
                self.upstream = int(a)
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
            elif a in ["-o", "-r", "-u", "-d", "-s", "-wm", "-wl", "-ws", "-w"]:
                prev = a
            elif a == "-f":
                self.toFasta = True
            elif self.genelist is None:
                ext = os.path.splitext(a)[1]
                if ext in [".gtf", ".GTF"]:
                    self.genelist = GTFParser(a)
                elif ext in [".gff", ".gff3", ".GFF", ".GFF3"]:
                    self.genelist = GFFParser(a)
                else:
                    sys.stderr.write("Error: gene database is not a GTF or GFF file.\n")
                    return False
            else:
                if a[0] == '@':
                    with open(a[1:], "r") as f:
                        self.genes += [ g.strip() for g in f.read().split("\n") ]
                else:
                    self.genes.append(a)
        if os.path.isfile(self.reference):
            self.seqman = SequenceManager(self.reference)
            if self.genelist and self.genes:
                self.genelist.load()
                return True
            
            else:
                return self.usage()
        else:
            sys.stderr.write("Error: please specify genome reference file with -r option.\n")
            return False

    def banner(self):
        sys.stdout.write("""\x1b[1;36m***************************************************
* fengc.py - design primers for FENGC experiments *
***************************************************\x1b[0m
""")

    def usage(self, what=None):
        sys.stdout.write("""
\x1b[1mUsage: fengc.py [options] genesdb genes...\x1b[0m

where `genesdb' is a gene database in GFT or GFF format, and `genes' is one or
more gene names, or (if preceded by @) a file containing gene names, one per line.

\x1b[1mGeneral options:\x1b[0m

  -r R | Specify location of genome reference file (required).
  -o O | Write output to file O (default: standard output).
  -f   | Write target sequences to FASTA files.

\x1b[1mDesign options:\x1b[0m

  -u U | Number of bp upstream of TSS (default: {}).
  -d D | Number of bp downstream of TSS (default: {}).
  -s S | Number of bp upstream/downstream of regions of interest (default: {}).

\x1b[1mWeight options (see -h weight for details):\x1b[0m

  -wm W | Set weight for MT penalty (default: {}).
  -wl W | Set weight for region length when larger than target (default: {}).
  -ws W | Set weight for region length when smaller than target (default: {}).
  -w  W | Set both -wl and -ws to W.

""".format(self.upstream, self.downstream, self.field, self.weightmt, self.weightlen1, self.weightlen2))
        return False

    def run(self):
        gnames = []
        gcoords = []
        self.regsize  = self.upstream + self.downstream # Size of target region

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

        self.sequences = self.seqman.getSequences(gcoords)

        sys.stderr.write("\n\x1b[1mFinding oligos:\x1b[0m\n")
        sys.stderr.write("  \x1b[4mGene\x1b[0m                \x1b[4mSize\x1b[0m    \x1b[4m Oligo 1 \x1b[0m   \x1b[4m Oligo 2 \x1b[0m   \x1b[4m Oligo 3 \x1b[0m\n")
        sys.stderr.write("                              MT    %GC   MT    %GC   MT    %GC\n")
        regstart = self.field                # Start of target region
        regend   = self.field + self.regsize # End of target region
        for i in range(len(gnames)):
            seq = self.sequences[i]
            if gcoords[i][3] == '-':
                seq.seq = revcomp(seq.seq)
            seq.name = gnames[i] + " " + seq.name[1:]
            if self.toFasta:
                seq.write(gnames[i] + ".fa")
            best = self.findOptimalOligos(seq, regstart, regend)
            seq.extra = best
            sys.stderr.write("  {:20}{}bp   {:.1f}  {}%   {:.1f}  {}%   {:.1f}  {}%\n".format(
                gnames[i], best[1].start - best[0].end, 
                best[0].mt, int(100*best[0].gcperc), 
                best[1].mt, int(100*best[1].gcperc), 
                best[2].mt, int(100*best[2].gcperc)))
        sys.stderr.write("\n\x1b[1mWriting output::\x1b[0m\n")

    def findOptimalOligos(self, seq, regstart, regend):
        positions = seq.findVT()
        pos1 = findCandidatePositions(positions, self.ncandout, regstart, -1) + findCandidatePositions(positions, self.ncandin, regstart, 1)
        pos1.sort(key=lambda p: abs(regstart - p))
        pos2 = findCandidatePositions(positions, self.ncandout, regend, 1) + findCandidatePositions(positions, self.ncandin, regend, -1)
        pos2.sort(key=lambda p: abs(regend - p))
        oligos1 = [ findOptimalMT(seq.seq, start, 1, self.minlength, self.maxlength) for start in pos1 ]
        oligos2 = [ findOptimalMT(seq.seq, start, 1, self.minlength, self.maxlength) for start in pos2 ]

        # Find best pair
        best = None
        bestscore = 1e6
        for a in oligos1:
            for b in oligos2:
                if a and b:
                    # See if we can find Oligo3
                    c = findOptimalMT(seq.seq, b.start, -1, self.minlength, self.maxlength)
                    if c:
                        score = self.score(a, b)
                        if score < bestscore:
                            bestscore = score
                            best = [a, b, c]
        return best

    def score(self, a, b):
        s1 = (a.mt - 65)**2
        s2 = (a.mt - 65)**2
        s3 = (a.mt - b.mt)**2
        size = b.start - a.end - self.regsize
        if size > 0:
            s4 = self.weightlen1 * size**2
        else:
            s4 = self.weightlen2 * size**2
        return self.weightmt*(s1 + s2 + s3) + s4

### Let's get things started

if __name__ == "__main__":
    args = sys.argv[1:]
    M = Main()
    M.banner()
    if M.parseArgs(args):
        M.run()
    else:
        M.usage(args)

### Test with:

## ./fengc.py -f -r /ufrc/data/reference/icbr/GRCh38/GRCh38.fa /ufrc/data/reference/icbr/GRCh38/Homo_sapiens.GRCh38.95.gtf APOE KDM6A TLR4