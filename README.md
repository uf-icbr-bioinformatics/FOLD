# FOLD


FOLD (FENGC OLigonucletide Designer) is a program to design sets of FENGC primers [insert ref here]
for one or more target genes. Primers are selected by optimizing a set of local and global constraints.
Local constraints apply to each individual target gene, while global constraints aim at ensuring 
that the resulting PCR fragments have similar sizes.

## Prerequisites
This program requires `bedtools` to be installed and in PATH. If the `primer3-py` package (available from
https://pypi.org/project/primer3-py/) is installed, the program will use it to perform Tm calculations,
otherwise it will fall back on an internal, less accurate method.

## Usage

```bash
fengc.py [options] reference genes...
```

where `reference` is a FASTA file containing the genome reference, and `genes` is one or
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

## Tm options
  
Option | Description
-------|------------
  -tl L   | Set minimum allowed Tm (default: 62).
  -th H   | Set maximum allowed Tm (default: 68).
  -ol L,H | Choose oligos in the size range L - H (default 8,30).
  -gc L,H | Set GC range to L - H (default: 40,65).
  -tm3    | Use primer3 method to compute Tm (default: builtin method).
  -pt T   | Limit temperature for hairpin/heterodimers (default: 50).
  -pd D   | Limit Ds for hairpin/heterodimers (default: -9).
  -mr R   | Maximum number of repeat bases (default: 6).

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

See the [Amplicon design](#amplicon-design) and [Amplicon size weighing](#amplicon-size-weighing-example) sections for more details.

## Weight options

Option | Description
-------|------------
  -wm W | Set weight for Tm penalty (default: 1.0).
  -wl W | Set weight for region length when larger than target (default: 1.0).
  -ws W | Set weight for region length when smaller than target (default: 1.0).
  -w  W | Set both -wl and -ws to W.

See the [Optimization weights](#optimization-weights) section for details.

## Optimization options:

Option | Description
-------|------------
  -n N | Perform N rounds of global optimization (default: 1000000).
  -l   | Local only - do not perform global optimization.

## Amplicon design

These are the steps used by FOLD to design an amplicon starting from the
TSS position of a gene. For each gene, three oligo primers need to be 
identified (called OligoA, OligoB, and OligoC respectively in the following).

* The `desired amplicon` by default extends from 400bp upstream of the TSS
  to 100bp downstream of the TSS. These values can be changed with the `-u`
  and `-d` options respectively.

* The two options `-mo` and `-mu` represent how much the amplicon size is allowed
  to grow or shrink, in percentage. The default values are 15 and 5, which means
  that the amplicon size can range from 475bp to 575bp.

* The program looks for all potential oligos around the left and right boundaries
  of the desired amplicon, and ranks them according to the procedure described in
  [Optimization weights](#optimization-weights). This is done separately for OligoA, OligoB, and OligoC, and the
  set with the best score is chosen as the optimal triple (this can be changed
  later if global optimization is enabled).

* The program forces OligoA to be upstream of the TSS and OligoB to be downstream
  of it. To remove this constraing (and allow amplicons that do not cover the 
  TSS) add the `-nt` option.

## Optimization weights

This program tests all possible pairs of candidate oligos at the 5' and 3'
end of the target sequence, assigning each pair a score, and selects the
pair with the optimal score. The score is the sum of five components:

  A. Squared difference between Tm for Oligo 1 and 65.

  B. Squared difference between Tm for Oligo 2 and 65.

  C. Squared difference between Tm for Oligo 3 and 65.

  D. Square of the difference between the highest and lowest Tms.

  E. Amplicon size factor.

The first four components are multiplied by weight `wm`. Component E can be computed 
in two different ways, depending on the value of the `-S` argument. If the value is 
"s" (the default), E is the predicted amplicon size (i.e. the distance between the 
start of Oligo1 and the start of Oligo2) multiplied by `wl` if it is larger than 
the desired amplicon size, or by `ws` if it is smaller. 

If -S is "t", component E is the square of the distance between the start of Oligo1 
and the desired amplicon start, multipled by `wl`, plus the square of the distance 
between the start of Oligo2 and the desired amplicon end, multiplied by `ws`. See
the [Amplicon size weighing](#amplicon-size-weighing-example) section for an example.

The values of the five components are added together to generate the final score.
Higher weights mean that the corresponding component has more impact on the choice
of the optimal pair, while a weight of 0 means that the corresponding component
has no impact. Examples:

* To select oligos based only on the size of the amplified region, ignoring the MT,
  use `-wm 0`.

* To avoid producing an amplified region shorter than the desired one, use `-ws 1000`
  (or any other very large number). This will assign an extremely high penalty
  to regions shorter than the desired length



## Amplicon size weighing example

In this example, the gene TSS is at position 1,000. The desired amplicon extends
from 400bp upstream of the TSS to 100b downstream, for a total size of 500bp.
Therefore, t1 is at position 600 and t2 at position 1100.

```
         600                                  1000      1100
         t1                                   TSS       t2
         |              400bp                 |  100bp  |
---------[####################################+#########]--------
       ^                                               ^
       A                                               B
```

Imagine that OligoA is at position 570, and OligoB at position 1090, with wl=2
and ws=3. Component E of the weight is computed in the following ways.

If -S is "s":

  Predicted amplicon size is 520, which is larger than the desired size (500),
  therefore E = 520*2 = 1040.

If -S is "t":

  The distance between A and t1 is 30, and between B and t2 is 10. Therefore
  E = (30^2) * 2 + (10^2) * 2 = 900 * 2 + 100 * 3 = 1800 + 300 = 2100.
