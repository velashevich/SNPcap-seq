# SNPcap-seq
This pipeline is used to calculate directional and absolute bias statistics for samples dual-sequenced with SNPcap-Seq and whole-genome shotgun sequencing (WGS).

First the genotype likelihoods have to be generated for each paired library (.bam formatted) independently using a Samtools genotype likelihood model (-GL 1) in ANGSD (Korneliussen et al., 2014). 
