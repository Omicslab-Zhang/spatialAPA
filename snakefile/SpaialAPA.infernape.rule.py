
import pandas as pd

# ================
# 读取配置文件
# ================
configfile: "config.infernape.yaml"

SIF_IMAGE   = config["sif_image"]
WORKDIR     = config["workdir"]
SAMPLE_FILE = config["sample_table"]

# 读取样本表
sample_table = pd.read_csv(SAMPLE_FILE, sep="\t", header=None, names=["GSE", "GSM"])
samples = list(sample_table.itertuples(index=False, name=None))

rule all:
    input:
        expand(f"{WORKDIR}/{{gse}}/{{gsm}}/infernape_output/done.txt",
               zip, gse=[s[0] for s in samples], gsm=[s[1] for s in samples])

rule run_infernape:
    input:
        bam = f"{WORKDIR}/{{gse}}/{{gsm}}/possorted_genome_bam.bam",
        whitelist = f"{WORKDIR}/{{gse}}/{{gsm}}/barcode.txt"
    output:
        done = f"{WORKDIR}/{{gse}}/{{gsm}}/infernape_output/done.txt"
    benchmark:
        f"{WORKDIR}/{{gse}}/{{gsm}}/benchres.txt"
    log:
        f"{WORKDIR}/{{gse}}/{{gsm}}/infernape.log"
    threads: 8
    shell:
        """
        singularity exec --pwd / {SIF_IMAGE} \
        bash -c '
            source /infernape/bin/activate && \
            Rscript /run_infernape.R \
              -b {input.bam} \
              -c {input.whitelist} \
              -w {WORKDIR}/{wildcards.gse}/{wildcards.gsm}/infernape_output \
              -g hg38 -t {threads}
        ' > {log} 2>&1 && \
        touch {output.done}
        """
