### Classes for fengc

import sys
import random
import subprocess as sp
from datetime import datetime

from fengcUtils import revcomp, findCandidatePositions, ngc, mt, nlower

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

    def __str__(self):
        return "#<Oligo {}-{}>".format(self.start, self.end)

    def __repr__(self):
        return "#<Oligo {}-{}>".format(self.start, self.end)

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

class SequenceManager(object):
    reference = None

    def __init__(self, reference):
        self.reference = reference

    def getSequences(self, genes):
        """Retrieve DNA sequence for each Sequence in the supplied genes."""
        wanted = {}
        proc = sp.Popen("bedtools getfasta -fi {} -bed /dev/stdin".format(self.reference), shell=True, stdin=sp.PIPE, stdout=sp.PIPE)

        # Pass coordinates in BED format to bedtools, and create dictionary to match
        # results with Sequence objects
        for g in genes:
            intv = g.transcripts[0]
            wanted[intv.label] = intv
            proc.stdin.write("{}\t{}\t{}\n".format(intv.chrom, intv.start, intv.end).encode('utf-8'))
        proc.stdin.close()

        while True:
            name = proc.stdout.readline()
            if not name:
                break
            name = name.decode('utf-8').rstrip("\n")[1:]
            seq = proc.stdout.readline().decode('utf-8').rstrip("\n")
            if name in wanted:
                s = wanted[name]
                if s.strand == "-":
                    seq = revcomp(seq)
                s.seq = seq
            else:
                sys.stderr.write("Unexpected output from bedtools: {}\n".format(name))

        good = []
        bad = []
        
        for g in genes:
            if g.transcripts[0].seq:
                good.append(g)
            else:
                bad.append(g)

        return good, bad

# Intervals

class Interval(object):
    name = ""                   # Source of this interval, e.g. a gene name
    label = ""                  # chr:start-end
    chrom = ""
    start = 0
    end = 0
    strand = "+"

    def __init__(self, name, chrom, start, end, strand="+"):
        self.name = name
        self.chrom = chrom
        self.start = start 
        self.end = end
        self.strand = strand
        self.label = "{}:{}-{}".format(chrom, start, end)

    def __str__(self):
        return "#<Interval {} {}>".format(self.name, self.label)

    def clone(self):
        return Interval(self.name, self.chrom, self.start, self.end, self.strand)

    def length(self):
        return self.end - self.start

# A Sequence is an Interval with associated sequence and other properties

class Sequence(Interval):
    seq = ""
    valid = True                # Set to False if no valid oligos found
    triples = []
    rndtriple = None            # One triple picked at random
    best = None                 # Best triple found so far
    bestGC = 0                  # GC content corresponding to the best triple

    def __init__(self, name, chrom, start, end, strand="+"):
        self.name = name
        self.chrom = chrom
        self.start = start 
        self.end = end
        self.strand = strand
        self.label = "{}:{}-{}".format(chrom, start, end)
        self.triples = []

    def __str__(self):
        return "#<Sequence {} {} {}bp>".format(self.name, self.label, len(self.seq))

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

    def findBestGC(self):
        if not self.bestGC:
            start = self.best.oligo1.start+1
            end = self.best.oligo2.start+1
            gc = ngc(self.seq[start:end])
            self.bestGC = 1.0 * gc / (end-start)
        return self.bestGC

    def findOptimalOligos(self, regstart, regend, main):
        tsspos = regstart + main.upstream
        positions = self.findVT()

        #sys.stderr.write("{}\n".format((regstart, regend, tsspos)))
        pos1 = findCandidatePositions(positions, main.ncandout, regstart, -1) + findCandidatePositions(positions, main.ncandin, regstart, 1)
        pos1.sort(key=lambda p: abs(regstart - p))
        #sys.stderr.write("{}\n".format(pos1))
        pos2 = findCandidatePositions(positions, main.ncandout, regend, 1) + findCandidatePositions(positions, main.ncandin, regend, -1)
        pos2.sort(key=lambda p: abs(regend - p))
        if main.includeTSS:
            pos1 = [ p for p in pos1 if p < tsspos ]
            pos2 = [ p for p in pos2 if p > tsspos ]

        #sys.stderr.write("{}+{} positions\n".format(len(pos1), len(pos2)))
        #sys.stderr.write("{}\n".format(pos2))
        oligos1 = [ self.findOptimalMT(start, 1, main) for start in pos1 ]
        oligos2 = [ self.findOptimalMT(start, 1, main) for start in pos2 ]

        for a in oligos1:
            for b in oligos2:
                if a and b:
                    # See if we can find Oligo3
                    c = self.findOptimalMT(b.start, -1, main)
                    if c:
                        t = Triple(a, b, c)
                        if main.regmin <= t.size <= main.regmax:
                            t.score = main.score(t)
                            self.triples.append(t)
        #sys.stderr.write("{} triples\n".format(len(self.triples)))
        self.triples.sort(key=lambda t: t.score)
        self.triples = self.triples[:main.maxtriples]

    # Find optimal oligo length in a sequence giving MT as close as possible to 65 degrees
    # Returns an Oligo object

    def findOptimalMT(self, start, direction, main):
        # print (start, direction, minlen, maxlen)
        bestlen = 0
        bestmt = 100
        delta = 100
        for l in range(main.minlength, main.maxlength):
            if direction == 1:
                oligo = self.seq[start:start+l] ## TODO: optimize - we just need to look at the base being added
            else:
                oligo = self.seq[start-l+1:start+1]
            gc = ngc(oligo)
            gcperc = 100.0 * gc / l

            # We want GC% between 40 and 65
            if gcperc < main.mingc or gcperc > main.maxgc:
                continue

            # Count lower-case (i.e. repeatmasked) bases
            nl = nlower(oligo)
            if nl >= main.maxrpt:
                continue

            # Check hairpin formation
            if main.mt_primer3:
                hairpin_info = main.calcHairpin(oligo)
                if hairpin_info.structure_found:
                    if hairpin_info.ds <= main.pr3_ds and hairpin_info.tm >= main.pr3_tm:
                        continue

            # Compute MT
            if main.mt_primer3:
                thismt = main.calcTm(oligo)
            else:
                thismt = mt(gc, l)

            # If we're closer to 65 than the current best, keep this
            if abs(thismt-65) < delta:
                bestlen = l
                bestmt = thismt
                delta = abs(thismt-65)

        if main.minmt <= bestmt <= main.maxmt:
            return Oligo(self.seq, start, direction, bestlen, bestmt)
        # else:
        return None

# Simple Gene DB

class Gene(object):
    accession = ""
    name = ""
    transcripts = []

    def __init__(self, accession, name):
        self.accession = accession
        self.name = name
        self.transcripts = []

    def __str__(self):
        return "#<Gene {}:{}>".format(self.name, self.accession)

    def clone(self):
        return Gene(self.accession, self.name)

    def addTranscript(self, accession, chrom, start, end, strand):
        """Add a new transcript with the given properties to this gene, unless 
one with the same TSS already exists."""
        for tr in self.transcripts:
            if (strand == '+' and start == tr.start) or (strand == '-' and end == tr.end):
                return 0          # already present, ignore
        self.transcripts.append(Interval(self.name + ":" + accession, chrom, start, end, strand=strand))
        return 1

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

    def badChroms(self, badchroms):
        self.out.write("*** Could not find reference sequence for the following {} genes:\n\n".format(len(badchroms)))
        for b in badchroms:
            self.out.write(b.name + "\n")
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

