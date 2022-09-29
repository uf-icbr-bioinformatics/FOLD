### Utilities for fengc

import subprocess as sp

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
    # else:
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

    i = 0
    # Find first position after `start'
    for i, p_i in enumerate(positions):
        if p_i >= start:
            break
    if direction == 1:
        return positions[i:i+howmany]
    # else:
    p = max(i-howmany, 0)
    return positions[p:i]
