import vcf
import numpy as np


def get_VAF_dict(records):
    cp_alleles_GT_VAF = dict()
    correct_haploid_GT_values = {'0/0', '0|0', '0/1', '0|1', '1/1', '1|1', './.', '0/.'}

    for row in records:
        chr_pos = f'{row.CHROM} {row.POS}'
        alleles = ','.join([str(x) for x in row.alleles])
        if chr_pos not in cp_alleles_GT_VAF:
            cp_alleles_GT_VAF[chr_pos] = dict()

        samples = row.samples

        for sample in samples:
            s_name = sample.sample
            AD = sample.data.AD
            # don't use the sample.data.DP as sum of AD is not equal to DP
            # -> it means that DP contains information by all reads (not only informative as in AD)
            DP = sum(AD)
            GT = sample.data.GT

            if GT in correct_haploid_GT_values and DP > 0:
                if alleles not in cp_alleles_GT_VAF[chr_pos]:
                    cp_alleles_GT_VAF[chr_pos][alleles] = dict()

                if GT not in cp_alleles_GT_VAF[chr_pos][alleles]:
                    cp_alleles_GT_VAF[chr_pos][alleles][GT] = list()

                cp_alleles_GT_VAF[chr_pos][alleles][GT].append([x / DP for x in AD])

    return cp_alleles_GT_VAF


def print_output(cp_alleles_GT_VAF):
    for chr_pos, alleles_GT_VAF in cp_alleles_GT_VAF.items():
        for alleles, GT_VAF in alleles_GT_VAF.items():
            for GT, VAF in GT_VAF.items():
                mx = np.matrix(VAF)
                lst_sum = list(np.array(mx.sum(axis=0))[0])
                lst_VAF = [x / mx.shape[0] for x in lst_sum]
                print(chr_pos, alleles, ','.join([str(x) for x in lst_VAF]))


if __name__ == "__main__":
    vcf_path = './vcf_data/Pf_M76611.pf7.200_samples.vcf'
    vcf_records = vcf.Reader(open(vcf_path, 'r'))

    cp_alleles_GT_VAF = get_VAF_dict(records=vcf_records)

    print_output(cp_alleles_GT_VAF=cp_alleles_GT_VAF)
