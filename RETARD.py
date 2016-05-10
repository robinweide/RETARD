# Make 5 whole genome reference-fastas (RWalpha, RWbeta, RWgamma, RWdelta, RWepsilon)
## Get random 500bp, 370bp and 156bp chunks of the whole genome
## Shuffle all chunks, so that the don't overlap
## Get GC-content of every chunk

# Make 5 acrocentric reference-fastas (RAalpha, RAbeta, RAgamma, RAdelta, RAepsilon)
## Get random 500bp, 370bp and 156bp chunks of the acrocentric chromosomes
## Shuffle all chunks, so that the don't overlap
## Get GC-content of every chunk

# Create shuffled rDNA-fastas (Ralpha, Rbeta, Rgamma, Rdelta, Repsilon)
## Chunk into 500bp, 370bp and 156bp
## Shuffle the chinks within the subunits
## Get GC-content of every chunk

# Map
## Extract fastq from input-bam
## Map to references

# Analysis part 1
## Get sum, mean, SD, median and MAD of the coverage per chunk (bedtools coverage -d)
## Output table: chunk | ref | a/b/g/d/e | mean GC | sumCoverage | meanCoverage | SDCoverage | medianCoverage | MADCoverage

# Analysis part 2
## Perform LOESS smoothed repeated running quantile on both types of reference-counts
## Fit models on rDNA-counts
## Output-table: rNDA-chunk | raw count | gc-corrected count (WG) | gc-corrected count (acrocentric)
