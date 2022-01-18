#!/bin/bash
cat - \
  | awk '$1~/chr[0-9X,^_]+$/' \
  | sed 's/chr//'
