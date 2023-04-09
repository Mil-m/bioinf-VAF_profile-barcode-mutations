import vcf
import numpy as np


def get_vaf_dict(records: vcf.Reader) -> dict:
    cp_alleles_gt_vaf = dict()
    correct_haploid_gt_values = {'0/0', '0|0', '0/1', '0|1', '1/1', '1|1', './.', '0/.'}

    for row in records:
        chr_pos = f'{row.CHROM} {row.POS}'
        alleles = ','.join([str(x) for x in row.alleles])
        if chr_pos not in cp_alleles_gt_vaf:
            cp_alleles_gt_vaf[chr_pos] = dict()

        samples = row.samples

        for sample in samples:
            ad_values = sample.data.AD
            # don't use the sample.data.DP as sum of AD is not equal to DP
            # -> it means that DP contains information by all reads (not only informative as in AD)
            dp_value = sum(ad_values)
            gt_value = sample.data.GT

            if gt_value in correct_haploid_gt_values and dp_value > 0:
                if alleles not in cp_alleles_gt_vaf[chr_pos]:
                    cp_alleles_gt_vaf[chr_pos][alleles] = dict()

                if gt_value not in cp_alleles_gt_vaf[chr_pos][alleles]:
                    cp_alleles_gt_vaf[chr_pos][alleles][gt_value] = list()

                cp_alleles_gt_vaf[chr_pos][alleles][gt_value].append([x / dp_value for x in ad_values])

    return cp_alleles_gt_vaf


def print_output(cp_alleles: dict, outp_file: str) -> None:
    outp_data = ''
    for chr_pos, alleles_gt_vaf in cp_alleles.items():
        for alleles, gt_vaf in alleles_gt_vaf.items():
            for gt, vaf in gt_vaf.items():
                mx = np.matrix(vaf)
                lst_sum = list(np.array(mx.sum(axis=0))[0])
                lst_vaf = [x / mx.shape[0] for x in lst_sum]

                data_vaf = ','.join([str(x) for x in lst_vaf])
                outp_data += f"{chr_pos} {alleles} {data_vaf}\n"

    with open(outp_file, 'w+') as f:
        f.write(outp_data)


if __name__ == "__main__":
    vcf_path = './vcf_vaf_data/Pf_M76611.pf7.200_samples.vcf'
    vcf_records = vcf.Reader(open(vcf_path, 'r'))

    cp_alleles_dict = get_vaf_dict(records=vcf_records)

    print_output(cp_alleles=cp_alleles_dict, outp_file='./vcf_vaf_data/outp_vaf.txt')
