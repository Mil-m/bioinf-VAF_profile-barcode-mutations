# bioinf-VAF_profile-barcode-mutations
Variant allele frequencies calculation and determination of a profile/barcode for mutations

`pip3 install -r requirements.txt`<br>

# Part A

input/output directory -> vcf_vaf_data directory<br>

(A1) output file there -> vcf_vaf_data/outp_vaf.txt<br>
these output data have duplicates because they have different genotype values for this haploid genome, for example<br>
`Pf_M76611 685 TTCGTACA,CTCGTACA,*,T,ATCGTACA,GTCGTACA 1.0,0.0,0.0,0.0,0.0,0.0 -- 0/0`<br>
`Pf_M76611 685 TTCGTACA,CTCGTACA,*,T,ATCGTACA,GTCGTACA 0.06818181818181818,0.9318181818181818,0.0,0.0,0.0,0.0 -- 1/1`<br>
If these differences can be explained by different sample conditions or other factors not related to mutational processes, then these data might be merged

python script -> calculate_VAF.py

# Part B

input directory -> mutations_data

ipython notebook -> mutations.ipynb

(B1) Output data in this notebook<br>
(B2) Answer to the question: The sample from “unknown_sample_mutations.txt” is more likely to have been sampled in August

