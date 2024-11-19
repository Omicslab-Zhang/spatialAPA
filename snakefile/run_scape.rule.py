

# Snakefile

rule scape_apa:
    input:
        bam = config['IN_PATH'] + "/{GSE}/{GSM}/possorted_genome_bam.bam",
        cb = config['IN_PATH'] + "/{GSE}/{GSM}/barcodes.tsv"
    output:
        pasite = config['IN_PATH'] + "/{GSE}/{GSM}/pasite.csv.gz"
    threads: 16
    params:
        outdir = config['IN_PATH'] + "/{GSE}/{GSM}/",
        cores = 16,
        bed = bed_file,
        scape_path = scape_path
    log:
        config['IN_PATH'] + "/{GSE}/{GSM}_scape.log"
    shell:
        """
        python {params.scape_path}/main.py apamix \
            --bed {params.bed} \
            --bam {input.bam} \
            --out {params.outdir} \
            --cores {params.cores} \
            --cb {input.cb} >{log} 2>&1
        """


rule step2:
    input:
        pasite = config['IN_PATH'] + "/{GSE}/{GSM}/pasite.csv.gz"
    output:
        collapse_pa = config['IN_PATH'] + "/{GSE}/{GSM}/collapse_pa.tsv.gz"
    params:
        labels="{GSM}",
        scape_path = config['scape_path']
    log:
        config['IN_PATH'] + "/{GSE}/{GSM}_scape_step2.log"
    shell:
        """
        python {params.scape_path}/scripts/group_pa.py \
            --files {input.pasite} \
            --labels {params.labels} \
            --outfile {output.collapse_pa} >{log} 2>&1
        """


rule step3:
    input:
        collapse_pa = config['IN_PATH'] + "/{GSE}/{GSM}/collapse_pa.tsv.gz",
    output:
        tmp_file_1 = config['IN_PATH'] + "/{GSE}/{GSM}/finish_step3.txt",
    params:
        scape_path=config['IN_PATH'],
        gse="{GSE}",
        gsm="{GSM}",
        scape_preprocess_1 = config['scape_preprocess_1'],
        gtf = config['gtf_file'],
        ensg_symbol = config['ensg_symbol']
    log:
        config['IN_PATH'] + "/{GSE}/{GSM}_scape_step3.log"
    shell:
        """
        Rscript {params.scape_preprocess_1} {params.scape_path} {params.gse} {params.gsm} {params.gtf} {params.ensg_symbol} \
            > {log} 2>&1 && touch {output.tmp_file_1}
        """


rule apalength:
    input:
        tmp_file_1 = config['IN_PATH'] + "/{GSE}/{GSM}/finish_step3.txt",
    output:
        length_csv = config['IN_PATH'] + "/{GSE}/{GSM}/gene_apa_length_mtx.csv"
    params:
        labels="{GSM}",
        inputdir = config['IN_PATH'],
        sample_list = config['gse_gsm_list'],
        stopcodon = config['stop_codon_file'],
        scape_calculate_len = config['calculate_apalen']
    log:
        config['IN_PATH'] + "/{GSE}/{GSM}_scape_step4.log"
    shell:
        """
        python {params.scape_calculate_len} \
            -i {params.inputdir} \
            -s {params.sample_list} \
            -n {params.labels} \
            -sc {params.stopcodon} \
            -o {output.length_csv} >{log} 2>&1
        """


rule apa_integrate:
    input:
        length_csv = config['IN_PATH'] + "/{GSE}/{GSM}/gene_apa_length_mtx.csv"
    output:
        tmp_file_2 = config['IN_PATH'] + "/{GSE}/{GSM}/finish_final.txt",
    params:
        scape_path=config['IN_PATH'],
        gse="{GSE}",
        gsm="{GSM}",
        scape_preprocess_2 = config['scape_preprocess_2']
    log:
        config['IN_PATH'] + "/{GSE}/{GSM}_scape_step5.log"
    shell:
        """
        Rscript {params.scape_preprocess_2} {params.scape_path} {params.gse} {params.gsm} > {log} 2>&1 && touch {output.tmp_file_2}
        """

