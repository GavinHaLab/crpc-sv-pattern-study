# Procedures
Use the following command for SRB detection. Bin size was set to 100kb. SV
events located within centromeric regions were excluded.
```
Rscript fishhook.r \
  -c \
  --no-centromere \
  -b 100 \
  -f none \
  --out-list srb.result.txt \
  --out-plot srb.qqplot.pdf \
  -e 0.75 \
  -i sv_data.bedpe
```
