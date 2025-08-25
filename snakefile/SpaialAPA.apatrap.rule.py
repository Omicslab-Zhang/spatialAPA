import pandas as pd

# ================
# 读取配置
# ================
configfile: "config.apatrap.yaml"

WORKDIR     = config["workdir"]
SAMPLE_FILE = config["sample_table"]
samtools_path = config["samtools"]
umitools_path = config["umitools"]
featureCounts_path = config["featureCounts"]
star_path = config["star"]
apatrap_R = config['spatialAPA_path']

# 读取样本表
sample_table = pd.read_csv(SAMPLE_FILE, sep="\t", header=None, names=["GSE", "GSM"])
samples = list(sample_table.itertuples(index=False, name=None))

rule all:
    input:
        expand(f"{WORKDIR}/{{gse}}/{{gsm}}/peak_counts.txt",
               zip, gse=[s[0] for s in samples], gsm=[s[1] for s in samples])

rule run_apatrap:
    input:
        bam = f"{WORKDIR}/{{gse}}/{{gsm}}/possorted_genome_bam.bam",
        whitelist = f"{WORKDIR}/{{gse}}/{{gsm}}/barcode.txt"
    output:
        done = f"{WORKDIR}/{{gse}}/{{gsm}}/done.txt",
        out_peak = f"{WORKDIR}/{{gse}}/{{gsm}}/peak_counts.txt"
    benchmark:
        f"{WORKDIR}/{{gse}}/{{gsm}}/benchres.txt"
    params:
        in_out_dir = f"{WORKDIR}/{{gse}}/{{gsm}}/",
    log:
        f"{WORKDIR}/{{gse}}/{{gsm}}/scAPAtrap.log"
    threads: 8
    shell:
        """
        Rscript {apatrap_R}/scripts/scapatrap.r {params.in_out_dir} {samtools_path} {umitools_path} {featureCounts_path} {star_path} > {log} 2>&1 && \
        touch {output.done}
        """
