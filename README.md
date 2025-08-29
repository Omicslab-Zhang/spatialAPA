# What is the spatial APA data analysis pipeline?
To quantify and integrate spatial APA data with gene expression, we developed a comprehensive spatial APA data analysis pipeline, called spatialAPA, built upon [Space Ranger](https://www.10xgenomics.com/cn/support/software/space-ranger/latest), [SCAPE](https://github.com/LuChenLab/SCAPE), and [Snakemake](https://snakemake.readthedocs.io/en/v7.32.4/). This pipeline encompasses multiple steps, including:
(1) Quality control of raw sequencing data (FASTQ/FQ files with corresponding tissue slice images).
(2) Quantification of APA events and gene expression.
(3) Calculation of Polyadenylation Site Index (PSI) and Weighted 3' UTR Length (WUL) scores for APA events.
(4) Downstream analysis: including cluster, spatial variable feature analysis, and GO/KEGG enrichment analysis.
(5) APA sites correction for single or cross sample analysis.

Users can easily run spatialAPA after configuring the required environment for their raw data ( fastq/fq of 10x visium spatial transcripomics data with attachment slice). The outputs, including APA matrices, PSI, WUL, and gene expression data, are seamlessly integrated into a Seurat object in R. This integration enables downstream spatial transcriptomics analyses using Seurat's extensive toolkit.

Using the spatialAPA pipeline, we constructed a comprehensive spatially resolved atlas of APA, analyzing 804,276 transcriptomic spots from 363 sections, spanning 56 projects, 18 organs, and 76 disease states. The resulting atlas provides unprecedented insights into the spatial dynamics of APA and is available for exploration at [SpatialAPAdb](http://www.biomedical-web.com/spatialAPAdb/home).

<div align="center">
  <img src="https://github.com/Omicslab-Zhang/spatialAPA/blob/main/image/figure1.tif" alt="Figure 1. Workflow of spatialAPAdb">
  <p><strong>Figure 1. Workflow of spatialAPA.</strong></p>
</div>

## Usage of the spatialAPA pipeline
The **[spaceranger](https://www.10xgenomics.com/cn/support/software/space-ranger/latest)** and **[SCAPE](https://github.com/LuChenLab/SCAPE)** tool should be pre-installed before using this pipeline. Other dependancy of environment for running pipeline of spatialAPA are list in environment.yml file. To ensure successful execution of the pipeline, **all configurations and essential parameters should be properly defined in the config.yaml file under snakefile directory.** You can download the [fastqs](https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR173/084/SRR17375084/) and [image data](https://github.com/Omicslab-Zhang/spatialAPA/tree/main/image/tissue_hires_image.png) for the example sample by clicking the link provided.
```
git clone https://github.com/Omicslab-Zhang/spatialAPA.git

cd spatialAPA

conda env create -f environment.yml

snakemake -s /path/to/spatialAPA/snakefile/spatialAPA.pipeline.py --configfile /path/to/spatialAPA/snakefile/config.yaml \
  -j 80 -k --scheduler greedy --rerun-incomplete > /path/to/log/file/spatialAPA.log 2>&1
```

## Explanation of Key Command-Line Options for Running the Spatial APA Pipeline
- **-s /path/to/spatialAPA/snakefile/spatialAPA.pipeline.py:**  
Specifies the Snakemake workflow file (Snakefile) containing the rules and steps for the pipeline. Replace /path/to/spatialAPA/snakefile/spatialAPA.pipeline.py with the actual file path of your pipeline.
- **--configfile /path/to/spatialAPA/snakefile/config.yaml**  
Specifies the YAML configuration file for the pipeline. This file contains all the necessary parameters and settings, such as input and output directories, reference files, computational resources, and other user-defined options.
By centralizing all these parameters in the config.yaml file, the workflow becomes modular and customizable.
- **-j 80:**  
Sets the maximum number of parallel jobs to run. In this example, the pipeline uses up to 80 threads to optimize computational efficiency.
- **-k:**  
Enables the pipeline to continue running independent jobs even if one job fails, ensuring partial progress in case of errors.
- **--scheduler greedy:**  
Utilizes the greedy scheduling algorithm to maximize resource utilization during the workflow execution.
- **--rerun-incomplete:**  
Ensures that any interrupted or partially completed jobs from previous runs are rerun to achieve a fully complete workflow.
- **/path/to/log/file/spatialAPA.log 2>&1:**  
Redirects both standard output and error messages from the pipeline execution to a specified log file for debugging and progress tracking. Replace /path/to/log/file/spatialAPA.log with the desired file path for your log file.
For additional options and detailed instructions, refer to the [Snakemake documentation](https://snakemake.readthedocs.io/en/v7.32.4/) or the pipeline's user guide.

## Merge polyA sites by windows size (default:±25nt)
This script is used to merge polyA sites from spatial transcriptomics or single-cell APA expression matrices. It helps to normalize genomic coordinates of APA events across different samples or tools by collapsing sites within a ±N nt window (default: ±25 nt). Each APA site should be named using the format: TSPAN8-chr12:71125093:-
```
python /path/to/script/spatialAPA_merge_pas.py \
    --sample_info /path/to/sample_info.txt \
    --window 25 \
    --outdir /path/to/out_directory
```
- **--sample_info: Path to the sample info file (required)**
- **--window: Window size for merging APA sites (default: 25)**
- **--outdir: Output directory to save merged APA matrix**

## Other options of spatialAPA pipeline
```
snakemake -s /path/to/spatialAPA/snakefile/SpaialAPA.apatrap.rule.py --configfile /path/to/spatialAPA/snakefile/config.apatrap.yaml \
  -j 80 -k --scheduler greedy --rerun-incomplete > /path/to/log/file/spatialAPA_apatrap.log 2>&1

snakemake -s /path/to/spatialAPA/snakefile/SpaialAPA.infernape.rule.py --configfile /path/to/spatialAPA/snakefile/config.infernape.yaml \
  -j 80 -k --scheduler greedy --rerun-incomplete > /path/to/log/file/spatialAPA_infernape.log 2>&1
```

## Further Information
To customize or debug the pipeline:
- Modify config.yaml to set appropriate paths, parameters, and sample-specific details.
- Check Snakemake's official documentation for advanced options and troubleshooting.
- Ensure necessary dependencies and paths are set correctly, as defined in the spatialAPA.pipeline.py.
- For more help, inspect specific rules or contact the pipeline maintainer.
