#!/bin/env python
"""Extract corresponding bands from a band file based on input table of bands
of interest. All sub-bands will be merged.
Input 1: cytoband files
<chr> <start> <end> <band_id>

Input 2: bands of interest table
<site_name> <band_id>
"""
import sys

class GenomicRegion:
  def __init__(self, chr, start, end, *args):
    self.chr = chr
    self.start = start
    self.end = end
    self.data = args

  def __str__(self):
    return "\t".join((self.chr, str(self.start), str(self.end), \
        "\t".join(self.data)))

  @classmethod
  def merge(cls, name, band, *args):
    """Merge intervals. Use very simple process, assuming all input regions
    *should* be merged. Basically just take the min start and max end.
    """
    chr = args[0].chr
    start = min([r.start for r in args])
    end = max([r.end for r in args])
    return cls(chr, start, end, name, band)

  @classmethod
  def from_str(cls, line):
    f = line.split()
    chr, start, end = f[0], int(f[1]), int(f[2])
    return cls(chr, start, end, *f[3:])

def main():
  if len(sys.argv) < 3:
    sys.stdout.write("Usage: %s <band_file> <table>\n"%sys.argv[0])
    sys.exit(1)

  band_pool = {}
  with open(sys.argv[1], 'r') as band_file:
    for l in band_file:
      f = l.split()
      band_pool[f[0][3:]+f[-1]] = GenomicRegion.from_str(l)

  with open(sys.argv[2], 'r') as list_file:
    for l in list_file:
      name, band = l.split()
      band_list = [band_pool[b] for b in band_pool if b.startswith(band)]
      print(GenomicRegion.merge(name, band, *band_list))

if __name__ == "__main__":
  main()
