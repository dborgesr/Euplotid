

```python
#Pipeline to go from ChiA-PET fastq reads all the way to interactions
!ChIA-PET2 -g /data/hg19.fa -f /data/fastq/SRR2054934_pass_1.fastq.gz -r /data/fastq/SRR2054934_pass_2.fastq.gz \
    -b /data/hg19.chrom.sizes -A ACGCGATATCTTATCTGACT -B AGTCAGATAAGATATCGCGT \
    -m 1 -s 1 -t 4 -M "-q .1" -k 2  -C 1
```
