#!/bin/env python
"""Convert the masked FASTA file to BED intervals.
Input:
>chr1
0000000000000000111111111111111111222222222111111

Output:
chr1	0	16	REG1	0
chr1	16	34	REG2	1
chr1	34	43	REG3	2
chr1	43	50	REG4	1
"""
import sys
from collections import defaultdict

def main():
  if len(sys.argv) < 2:
    sys.stdout.write("Usage: %s <masked_fa> > <mask_bed>\n"%sys.argv[0])
    sys.exit(1)

  chr_pool = defaultdict(list)
  with open(sys.argv[1], 'r') as inf:
    for l in inf:
      if l.startswith(">"):
        chr_name = l[1:].split()[0]
      else:
        chr_pool[chr_name].append(l.strip())

  for chr in chr_pool:
    seq = "".join(chr_pool[chr])
    start = 0
    for i in range(1, len(seq)):
      curr = seq[i]
      prev = seq[i-1]
      if curr != prev:
        sys.stdout.write("\t".join((chr, str(start), str(i), "mask", prev)) \
            + "\n")
        start = i
    sys.stdout.write("\t".join((chr, str(start), str(len(seq)), "mask", curr)) \
        + "\n")

if __name__ == "__main__":
  main()
