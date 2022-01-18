#!/bin/env python
"""Removing SV/CN data from input blacklist. Support * to exclude
an entire sample.
Input example:
  sample	chr
  DTB-205-BL	chr1
  DTB-023-BL	*
  ...
"""
import sys,re

def parse_list(inf):
  bl = {}
  for l in inf:
    if l.startswith("#"):
      pass
    f = l.split()
    sample, chr = f[0], f[1]
    if sample in bl:
      bl[sample].append(chr)
    else:
      bl[sample] = [chr]

  return bl

def main():
  if len(sys.argv) < 2:
    sys.stdout.write("Usage: %s <blacklist.txt>\n"%sys.argv[0])
    sys.exit()

  with open(sys.argv[1], 'r') as inf:
    bl = parse_list(inf)

  for l in sys.stdin:
    sample_check = [s if re.search(s, l) else False for s in bl]
    skip = False
    if any(sample_check):
      for s in filter(None, sample_check):
        chr_check = [c if c == "*" or re.search("\\b%s\\b"%c, l) else False \
            for c in bl[s]]
      if any(chr_check):
        skip = True
    if not skip:
      sys.stdout.write(l)

if __name__ == "__main__":
  main()
