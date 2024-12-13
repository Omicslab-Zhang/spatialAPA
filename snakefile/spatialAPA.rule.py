
# Snakefile
rule spaceranger_count:
    input:
        fastqs=config['IN_PATH'] + "/{GSE}/{GSM}/",
        image=lambda wildcards: wildcards.get("image", None),
        darkimage=lambda wildcards: wildcards.get("darkimage", None),
        loupe=lambda wildcards: wildcards.get("loupe", None),
    output:
        bam = config['IN_PATH'] + "/{GSE}/{GSM}/possorted_genome_bam.bam"
    params:
        transcriptome=config['REF_APP'],
        slide=config['SLIDE'],
        localcores=config['CORE'],
        localmem=config['MEM'],
        spatialAPA=config['spatialAPA_path'],
        outdir=config['IN_PATH'] + "/{GSE}"
    shell:
        """
        {params.spatialAPA}/scripts/spatialAPA_Gene.sh \
            {wildcards.GSM} \
            {input.fastqs} \
            {params.transcriptome} \
            {params.localcores} \
            {params.localmem} \
            {params.slide} \
            {params.outdir} \
            {input.image or ""} \
            {input.darkimage or ""} \
            {input.loupe or ""}
        """

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
        bed = config['bed_file'],
        scape_path = config['scape_path']
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

rule step3:
    input:
        pasite = config['IN_PATH'] + "/{GSE}/{GSM}/pasite.csv.gz"
    output:
        collapse_pa = config['IN_PATH'] + "/{GSE}/{GSM}/collapse_pa.tsv.gz"
    params:
        labels="{GSM}",
        scape_path = config['scape_path']
    log:
        config['IN_PATH'] + "/{GSE}/{GSM}_step3.log"
    shell:
        """
        python {params.scape_path}/scripts/group_pa.py \
            --files {input.pasite} \
            --labels {params.labels} \
            --outfile {output.collapse_pa} >{log} 2>&1
        """

rule step4:
    input:
        collapse_pa = config['IN_PATH'] + "/{GSE}/{GSM}/collapse_pa.tsv.gz",
    output:
        tmp_file_1 = config['IN_PATH'] + "/{GSE}/{GSM}/finish_step4.txt",
    params:
        res_path=config['IN_PATH'],
        gse="{GSE}",
        gsm="{GSM}",
        spatialAPA=config['spatialAPA_path'],
        gtf = config['gtf_file'],
        ensg_symbol = config['ensg_symbol']
    log:
        config['IN_PATH'] + "/{GSE}/{GSM}_step4.log"
    shell:
        """
        Rscript {params.spatialAPA}/scripts/spatialAPA_Annot.r \
            {params.res_path} {params.gse} {params.gsm} {params.gtf} {params.ensg_symbol} \
            > {log} 2>&1 && touch {output.tmp_file_1}
        """

rule apalength:
    input:
        tmp_file_1 = config['IN_PATH'] + "/{GSE}/{GSM}/finish_step4.txt",
    output:
        length_csv = config['IN_PATH'] + "/{GSE}/{GSM}/gene_apa_length_mtx.csv"
    params:
        labels="{GSM}",
        res_path = config['IN_PATH'],
        sample_list = config['gse_gsm_list'],
        stopcodon = config['stop_codon_file'],
        spatialAPA=config['spatialAPA_path'],
    log:
        config['IN_PATH'] + "/{GSE}/{GSM}_step5.log"
    shell:
        """
        python {params.spatialAPA}/scripts/spatialAPA_WUL.py \
            -i {params.res_path} \
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
        res_path=config['IN_PATH'],
        gse="{GSE}",
        gsm="{GSM}",
        spatialAPA=config['spatialAPA_path'],
    log:
        config['IN_PATH'] + "/{GSE}/{GSM}_integrate.log"
    shell:
        """
        Rscript {params.spatialAPA}/scripts/spatialAPA_APAPSI.r {params.res_path} {params.gse} {params.gsm} > {log} 2>&1 && touch {output.tmp_file_2}
        """

