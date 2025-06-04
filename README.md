# SNPcap-seq
## This pipeline is used for calculating and comparing directional and absolute bias statistics for samples dual-sequenced with SNPcap-Seq and whole-genome shotgun sequencing (WGS). The procedure relies on posterior genotype probability data.

First the genotype likelihoods have to be generated for each paired library (.bam formatted) independently using a Samtools genotype likelihood model (-GL 1) and only for the positions targeted by capture array (-rf, use a file chr:position format) using ANGSD (Korneliussen et al., 2014). The prior and posterior genotype probabilities can then be calculated following (Plassais et al., 2022).

Prior probabilities can be calculated using the simplified vcf file containing the coordinates of targeted SNP sites (should contain following columns: "Chromosome", "Position", "Ref", "Alt"). This has to be done only once for a specific capture panel. The script SS4.py is used to assign a probability of 0.31 to three genotypes compatible with the specific reference-alternative allele pair at a given site, while the other seven genotypes get a probability of 0.01:

python3 prior_generator_for_angsd_genolik.py <simplified_vcf_input> <output_name>

Posterior genotypes can then be computed with SS5.py using a Bayesian-based approach where we update the priors with linearized genotype likelihood ANGSD output:

python3 posterior_estimator_for_angsd_genolik.py <prior_file> <likelihood_file> <output_name>
