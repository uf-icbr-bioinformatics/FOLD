# FOLD

FOLD (FENGC OLigonucletide Designer) is a program to design sets of FENGC primers [insert ref here]
for one or more target genes. Primers are selected by optimizing a set of local and global constraints.
Local constraints apply to each individual target gene, while global constraints aim at ensuring 
that the resulting PCR fragments have similar sizes.

## Prerequisite
This program requires `bedtools` to be installed and in PATH. If the `primer3-py` package (available from
https://pypi.org/project/primer3-py/) is installed, the program will use it to perform Tm calculations,
otherwise it will fall back on an internal, less accurate method.

## Usage

```bash
fengc.py [options] reference genes...
```

where `reference' is a FASTA file containing the genome reference, and `genes' is one or
more gene names, or (if preceded by @) a file containing gene names, one per line.

## General options:

Option | Description
-------|------------
  -g G  | Use gene database G (in GTF or GFF format). Default: a file with the same name as the reference file, with extension .gtf or .gff (lower or uppercase).
  -c    | Add "chr" at the beginning of chromosome names that don't have it.
  -L    | For genes having multiple isoforms, select longest one (default: select all).
  -o O  | Write output to file O (default: standard output).
  -f    | Write target sequences to separate FASTA files.
  -F F  | Write target sequences to single FASTA file F.
  -O O  | Write oligo sequences to a FASTA file named O.
  -r R  | Write report to file R (default: report.txt).
  -D    | Check for potential heterodimers.

## Tm options:
  
Option | Description
-------|------------
  -tl L   | Set minimum allowed Tm (default: 62).
  -th H   | Set maximum allowed Tm (default: 68).
  -ol L,H | Choose oligos in the size range L - H (default 8,30).
  -gc L,H | Set GC range to L - H (default: 40,65).
  -tm3    | Use primer3 method to compute Tm (default: builtin method).
  -pt T   | Limit temperature for hairpin/heterodimers (default: 50).
  -pd D   | Limit Ds for hairpin/heterodimers (default: -9).
  -mr R   | Maximum number of repeat bases (default: 6)."

## Design options

Option | Description
-------|------------
  -u U  | Number of bp upstream of TSS (default: 400).
  -d D  | Number of bp downstream of TSS (default: 100).
  -s S  | Number of bp upstream/downstream of regions of interest (default: 2000).
  -mo O | Maximum % increase of amplicon size (default: 15%).
  -mu U | Maximum % decrease of amplicon size (default: 5%).
  -S A  | Amplicon sizing method, one of `s' or `t' (default: s).
  -nt   | Do NOT force amplicon to contain TSS.

see -h design and -h size for details):

## Weight options

Option | Description
-------|------------
  -wm W | Set weight for Tm penalty (default: 1.0).
  -wl W | Set weight for region length when larger than target (default: 1.0).
  -ws W | Set weight for region length when smaller than target (default: 1.0).
  -w  W | Set both -wl and -ws to W.

 (see -h weight for details)

## Optimization options:

Option | Description
-------|------------
  -n N | Perform N rounds of global optimization (default: 1000000).
  -l   | Local only - do not perform global optimization.

