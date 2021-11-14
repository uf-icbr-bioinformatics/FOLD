### Parsers for fengc

import sys
import csv

from fengcClasses import Gene

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

