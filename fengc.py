#!/bin/bash

import sys
import csv
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

# Extracting sequences for regions

class Sequence(object):
    name = ""
    seq = ""

    def __init__(self, name, seq):
        self.name = name
        self.seq = seq

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
            name = proc.stdout.readline()
            seq = proc.stdout.readline()
            seqs.append(Sequence(name, seq))
        return seqs
