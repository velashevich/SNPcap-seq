import pandas as pd
import numpy as np
import math
import sys
# This function rescales ANGSD genotype likelihoods and transform them into posterior
# It takes two lists as input: ANGSD genotype likelihood ratios and the prior
def do_post(l,p):
	num_l=[math.pow(10,float(value)) for value in l]
	num=[prior*likelihood for prior,likelihood in zip(num_l,p)]
	den=sum(num)
	pr=[likelihood/den for likelihood in num]
	return(pr)


# Importing the set of informed prior, one variant per row.
# Probability values corresponds to the genotypes in alphabetic order:
# AA, AC, AG, AT, CC, CG, CT, GG, GT, TT
prior_file=sys.argv[1]
prior_df = pd.read_csv(prior_file, header=0, skip_blank_lines=True, sep="\t")
prior_df["Prior"]=prior_df.apply(lambda row: list(row[3:]) , axis=1)
subset_prior=prior_df[["Variant_ID","Ref","Alt","Prior"]]
# Importing likelihhod data frame from file
likelihood_file=sys.argv[2]
outfile=sys.argv[3]
try:
	likelihood_df = pd.read_csv(likelihood_file, compression='gzip', header=None, skip_blank_lines=True, sep="\t")
except pd.errors.EmptyDataError:
	likelihood_df = pd.DataFrame()
if (likelihood_df.size>0):
	likelihood_df.columns = ["Chromosome","Position","AA","AC","AG","AT","CC","CG","CT","GG","GT","TT"]
	likelihood_df["lik"]=likelihood_df.apply(lambda row: list(row[2:]) , axis=1)
	likelihood_df["Variant_ID"]=likelihood_df.apply(lambda row: row.Chromosome.lstrip("chr") + "_" + str(row.Position), axis=1)
	subset_likelihood=likelihood_df[["Variant_ID","lik"]]
	# Merging Dataframes and calculating posterior
	final_df=subset_prior.merge(subset_likelihood, how="inner", on="Variant_ID")
	posterior_df=final_df.apply(lambda row: do_post(row.lik,row.Prior) , result_type="expand", axis=1)
	posterior_df.columns = ["AA","AC","AG","AT","CC","CG","CT","GG","GT","TT"]
	posterior_output=final_df[["Variant_ID","Ref","Alt"]].join(posterior_df)
else:
	posterior_output=pd.DataFrame(0, index=np.arange(0), columns=["Variant_ID","Ref","Alt","AA","AC","AG","AT","CC","CG","CT","GG","GT","TT"])
posterior_output.to_csv(outfile, header=True, index=False, sep="\t")
