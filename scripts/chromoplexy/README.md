# Input data
  - SV calls `sv_data.bedpe`
  - Copy-number calls `cn_data.seg`
  - Telomere region annotation `telomere_q5m.bed`
  - Chromothripsis calls `chromothripsis.bed`

# Other files
  - Parameter file for ChainFinder `parameters.txt`
  - Script for filtering CN segments `prepare_seg_from_TITAN.r`
  - Samples to be excluded from this analysis `blacklist.txt`
  - Script for filtering SV and CN calls `filter_blacklist.py`

# Procedures
See `run.sh`.
