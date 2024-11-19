# #-*- coding:utf-8 -*-
import os
import sys
import yaml

# 从命令行参数获取GSE值
#gse = config["GSE"]  ## 1
#gsm = config["GSM"]  ## 2

gse_gsm_infile = config["gse_gsm_list"]
# /public/home/jiangzh/projects/sc_apa/gse_gsm_101.txt

new_dict = {}
with open(gse_gsm_infile, 'r') as fl:
  for lineL in fl:
    line = lineL.strip().split("\t")
    if line[0].startswith("GSE ID"):
      continue
    else:
      gse_tmp = line[0]
      gsm_tmp = line[1]
      if gse_tmp in new_dict:
          new_dict[gse_tmp].append(gsm_tmp)
      else:
          new_dict[gse_tmp] = [gsm_tmp]

IN_PATH = config["IN_PATH"]
THREADS = 16

include: config["scape_rule"]
# /public/home/jiangzh/projects/sc_apa/run_scape.rule.py

# 展平GSE和GSM的组合
flattened_gse_gsm = [(gse, gsm) for gse in new_dict.keys() for gsm in new_dict[gse]]

rule all:
    input:
        expand(IN_PATH + "/{GSE}/{GSM}/pasite.csv.gz", zip, GSE=[gse for gse, gsm in flattened_gse_gsm], GSM=[gsm for gse, gsm in flattened_gse_gsm]),
        expand(IN_PATH + "/{GSE}/{GSM}/collapse_pa.tsv.gz", zip, GSE=[gse for gse, gsm in flattened_gse_gsm], GSM=[gsm for gse, gsm in flattened_gse_gsm]),
        expand(IN_PATH + "/{GSE}/{GSM}/finish_step3.txt", zip, GSE=[gse for gse, gsm in flattened_gse_gsm], GSM=[gsm for gse, gsm in flattened_gse_gsm]),
        expand(IN_PATH + "/{GSE}/{GSM}/gene_apa_length_mtx.csv", zip, GSE=[gse for gse, gsm in flattened_gse_gsm], GSM=[gsm for gse, gsm in flattened_gse_gsm]),
        expand(IN_PATH + "/{GSE}/{GSM}/finish_final.txt", zip, GSE=[gse for gse, gsm in flattened_gse_gsm], GSM=[gsm for gse, gsm in flattened_gse_gsm]),
