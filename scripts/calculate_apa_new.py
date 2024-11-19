import pandas as pd
import pysam
import gzip
import os
# import sys
import re
import argparse

'''
usage: 
python calculate_apa_new.py -n GSM6177599 \
    -o /public/home/jiangzh/projects/sc_apa/database/SCAPE_res/GSE203612/GSM6177599/gene_apa_length_mtx.csv
'''

# 样本信息
# samp = sys.argv[1]
# samp = "SPMDBS00000112"

# 读取PolyA位点文件并提取条形码
def read_pasite_csv(file_path):
    with gzip.open(file_path, 'rt') as f:
        df = pd.read_csv(f, index_col=0)
    return df

def calculateWUL(IN_PATH, gse_gsm_list, samp, stop_codon_file, output_file):
    # 定义文件目录
    rds_path = IN_PATH + "/spatial_h5ad/apa/rds"
    # 创建样本字典
    dict_file = gse_gsm_list
    dict_samp = {}
    with open(dict_file) as fl:
        for line in fl:
            lineL = line.strip().split()
            # dataset_id = lineL[0]
            # samp_id = lineL[1]
            gse_id = lineL[0]
            gsm_id = lineL[1]
            dict_samp[gsm_id] = {
                'gse_id': gse_id,
                'gsm_id': gsm_id                
            }
            # dict_samp[samp_id] = {
            #     'dataset_id': dataset_id,
            #     'gse_id': gse_id,
            #     'gsm_id': gsm_id
            # }

    # 加载APA注释文件
    '''
    示例文件:
        seqnames	pa_loc	strand	annot.start	annot.end	annot.tx_name	annot.gene_id	annot.symbol	annot.entrez	Type	Rank	pa_site	apa_id
    MFSD14A-chr1:100083332:+	chr1	100083332	+	100082210	100083377	ENST00000370152.8	ENSG00000156875.14	MFSD14A	NA	3UTRs	3UTRs	MFSD14A-chr1:100083332:+	MFSD14A-chr1:100083332:+
    RTCA-chr1:100291778:+	chr1	100291778	+	100291505	100292769	ENST00000370128.9	ENSG00000137996.13	RTCA	NA	3UTRs	3UTRs	RTCA-chr1:100291778:+	RTCA-chr1:100291778:+
    VCAM1-chr1:100739038:+	chr1	100739038	+	100738284	100739045	ENST00000294728.7	ENSG00000162692.12	VCAM1	NA	3UTRs	3UTRs	VCAM1-chr1:100739038:+	VCAM1-chr1:100739038:+
    UBE4B-chr1:10181242:+	chr1	10181242	+	10180368	10181367	ENST00000253251.12	ENSG00000130939.20	UBE4B	NA	LastExon1Kb	ExonRank27	UBE4B-chr1:10181242:+	UBE4B-chr1:10181242:+
    AMY2B-chr1:103579533:+	chr1	103579533	+	103579501	103579534	ENST00000684275.1	ENSG00000240038.8	AMY2B	NA	3UTRs	3UTRs	AMY2B-chr1:103579533:+	AMY2B-chr1:103579533:+
    KIF1B-chr1:10381628:+	chr1	10381628	+	10381604	10382603	ENST00000676179.1	ENSG00000054523.20	KIF1B	NA	LastExon1Kb	ExonRank49	KIF1B-chr1:10381628:+	KIF1B-chr1:10381628:+
    PGD-chr1:10420140:+	chr1	10420140	+	10419750	10420511	ENST00000270776.13	ENSG00000142657.21	PGD	NA	3UTRs	3UTRs	PGD-chr1:10420140:+	PGD-chr1:10420140:+
    PGD-chr1:10420509:+	chr1	10420509	+	10419750	10420511	ENST00000270776.13	ENSG00000142657.21	PGD	NA	3UTRs	3UTRs	PGD-chr1:10420509:+	PGD-chr1:10420509:+
    '''
    APA_file = os.path.join(rds_path, dict_samp[samp]['gse_id'], dict_samp[samp]['gsm_id'], "APA_annotations.csv")
    apa_annotation = pd.read_csv(APA_file, index_col=0)
    APAs = apa_annotation.index.values.tolist()

    # 加载BAM文件
    bam_file = os.path.join(IN_PATH, dict_samp[samp]['gse_id'], dict_samp[samp]['gsm_id'], 'possorted_genome_bam.bam')
    bam = pysam.AlignmentFile(bam_file, 'rb')  # pysam

    # 加载终止子文件
    '''
    示例文件:
    chr1	70006	70008	+	OR4F5
    chr1	450740	450742	-	OR4F29
    chr1	685716	685718	-	OR4F16
    chr1	944151	944153	+	SAMD11
    chr1	944694	944696	-	NOC2L
    chr1	963384	963386	+	KLHL17
    chr1	965189	965191	+	KLHL17
    chr1	974573	974575	+	PLEKHN1
    chr1	976172	976174	-	PERM1
    chr1	999059	999061	-	HES4
    '''
    #stop_codon_file = "/public/home/jiangzh/projects/sc_apa/GRCh38.v44.stop_codon.tsv"
    dict_gene_stop_codon = {}
    with open(stop_codon_file) as file:
        for line in file:
            lineL = line.strip().split()
            chr_st = lineL[0]
            start_st = int(lineL[1])
            end_st = int(lineL[2])
            strand_st = lineL[3]
            gene_st = lineL[4]
            if gene_st not in dict_gene_stop_codon:
                dict_gene_stop_codon[gene_st] = []
            dict_gene_stop_codon[gene_st].append({
                'chromosome': chr_st,
                'start': start_st,
                'end': end_st,
                'strand': strand_st
            })

    pasite_file = os.path.join(IN_PATH, dict_samp[samp]['gse_id'], samp, 'pasite.csv.gz')
    df = read_pasite_csv(pasite_file)
    barcodes = df.columns.tolist()

    # 计算每个spot的reads数目
    reads_per_spot = {barcode: {} for barcode in barcodes}
    for APA in APAs:
        # match = re.match(r"(.+)-chr(\d+|X|Y):(\d+):([+-])", APA)
        match = re.match(r"(.+)-(chr[\d+|X|Y]+):(\d+):([+-])", APA)
        if match:
            gene_apa, chr_apa, start_apa, strand_apa = match.groups()
            start_apa = int(start_apa)
            end_apa = start_apa + 1  # 假设APA位点为单碱基位置
        else:
            continue  # 如果匹配失败，跳过这个APA
        # 计算到最近stop_codon的距离
        dist_st = float('inf')
        tmp_st_list = dict_gene_stop_codon.get(gene_apa, [])
        for each in tmp_st_list:
            tmp_dist_st = abs(int(each['start']) - start_apa)
            if tmp_dist_st < dist_st:
                dist_st = tmp_dist_st
        # 提取APA位点对应的reads
        try:
            for read in bam.fetch(chr_apa, start_apa, end_apa):
                if read.has_tag('CB'):
                    barcode = read.get_tag('CB')
                    if barcode in barcodes:
                        if gene_apa not in reads_per_spot[barcode]:
                            reads_per_spot[barcode][gene_apa] = []
                        reads_per_spot[barcode][gene_apa].append({
                            'dist_st': dist_st,
                            'reads_count': 1
                        })
        except ValueError as e:
            print(f"Error processing APA: {APA} - {e}")
            continue

    # 计算WUL
    WUL_dict = {}
    for barcode, genes in reads_per_spot.items():
        for gene_apa, values in genes.items():
            if gene_apa not in WUL_dict:
                WUL_dict[gene_apa] = {}
            if barcode not in WUL_dict[gene_apa]:
                WUL_dict[gene_apa][barcode] = []
            total_reads = sum(item['reads_count'] for item in values)
            WUL_numerator = sum(item['dist_st'] * item['reads_count'] for item in values)
            WUL_value = WUL_numerator / total_reads if total_reads > 0 else 0
            WUL_dict[gene_apa][barcode].append(WUL_value)

    # 汇总成gene*spot矩阵
    gene_spot_matrix = {}
    for gene, spots in WUL_dict.items():
        for spot, WUL_values in spots.items():
            if gene not in gene_spot_matrix:
                gene_spot_matrix[gene] = {}
            gene_spot_matrix[gene][spot] = round(sum(WUL_values),2)
    
    gene_spot_df = pd.DataFrame.from_dict(gene_spot_matrix, orient='index')
    gene_spot_df.fillna(0, inplace=True)
    gene_spot_df.to_csv(output_file, sep=',')

def main():
    parser = argparse.ArgumentParser(description="Calculating WUL for SCAPE results.")
    parser.add_argument("-i", "--inputdir", help="The input path.")
    parser.add_argument("-s", "--samplelist", help="The input sample list.|project samplename|")
    parser.add_argument("-n", "--name", help="The input sample name (GSM).")
    parser.add_argument("-sc", "--stopcodon", help="The input stop codon file.")
    parser.add_argument("-o", "--out", help="The output file path.")
    args = parser.parse_args()
    calculateWUL(args.inputdir, args.samplelist, args.name, args.stopcodon, args.out)

if __name__ == "__main__":
    main()

