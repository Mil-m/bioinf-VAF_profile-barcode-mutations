import vcf
import numpy as np


def get_all_gt(records: vcf.Reader) -> list:
    """
    Receiving all genotype values from vcf.Reader records
    :param records: vcf.Reader records
    :return: list with all genotype values
    """
    all_gt_values = list()

    for row in records:
        samples = row.samples
        for sample in samples:
            gt_value = sample.data.GT
            if gt_value not in all_gt_values:
                all_gt_values.append(gt_value)

    return all_gt_values


def get_vaf_dict(records: vcf.Reader, correct_haploid_gt_values: set) -> dict:
    """
    Receiving cp_alleles_gt_vaf dictionary in the format below where
    sample.data.AD1 / sum(sample.data.AD) is a VAF for each sample):
        { chr_pos: alleles (for correct haploid genotype values):
            gt_value: [sample.data.AD1 / sum(sample.data.AD), ..., sample.data.ADn / sum(sample.data.AD)]
        }
        example:
        { 'Pf_M76611 33': 'T,A': '0/0': [[1.0, 0.0], ..., [1.0, 0.0]] }
        { 'Pf_M76611 33': 'T,A': '0|1': [[0.9166666666666666, 0.08333333333333333]] }
    :param records: vcf.Reader
    :param correct_haploid_gt_values: set with the correct haploid genotype values
    :return: cp_alleles_gt_vaf dictionary
    """
    cp_alleles_gt_vaf = dict()

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
            gt_value = sample.data.GT.replace('0|0', '0/0').replace('0|1', '0/1').replace(
                '1|1', '1/1').replace('0/.', './.')

            if gt_value in correct_haploid_gt_values and dp_value > 0:
                if alleles not in cp_alleles_gt_vaf[chr_pos]:
                    cp_alleles_gt_vaf[chr_pos][alleles] = dict()

                if gt_value not in cp_alleles_gt_vaf[chr_pos][alleles]:
                    cp_alleles_gt_vaf[chr_pos][alleles][gt_value] = list()

                cp_alleles_gt_vaf[chr_pos][alleles][gt_value].append([x / dp_value for x in ad_values])

    return cp_alleles_gt_vaf


def write_output(cp_alleles: dict, outp_file: str) -> None:
    """
    Receiving mean VAF for each sample set for genotype and writing in the output file
    :param cp_alleles: cp_alleles dictionary in the format below where
    sample.data.AD1 / sum(sample.data.AD) is a VAF for each sample):
        { chr_pos: alleles (for correct haploid genotype values):
            gt_value: [sample.data.AD1 / sum(sample.data.AD), ..., sample.data.ADn / sum(sample.data.AD)]
        }
    :param outp_file: output filepath
    :return: None
    """
    outp_data = ''
    for chr_pos, alleles_gt_vaf in cp_alleles.items():
        for alleles, gt_vaf in alleles_gt_vaf.items():
            for gt, vaf in gt_vaf.items():
                mx = np.matrix(vaf)
                sum_array = np.array(mx.sum(axis=0))[0]
                vaf_lst = [x / mx.shape[0] for x in sum_array]

                data_vaf = ','.join([str(x) for x in vaf_lst])
                outp_data += f"{chr_pos} {alleles} {data_vaf}\n"

    with open(outp_file, 'w+') as f:
        f.write(outp_data)


if __name__ == "__main__":
    vcf_path = './vcf_vaf_data/Pf_M76611.pf7.200_samples.vcf'

    all_gt_values = get_all_gt(records=vcf.Reader(open(vcf_path, 'r')))
    print(f"All genotype values list for the input records: {all_gt_values}")

    cp_alleles_dict = get_vaf_dict(
        records=vcf.Reader(open(vcf_path, 'r')),
        correct_haploid_gt_values={'0/0', '0|0', '0/1', '0|1', '1/1', '1|1', './.', '0/.'}
    )
    write_output(cp_alleles=cp_alleles_dict, outp_file='./vcf_vaf_data/outp_vaf.txt')
