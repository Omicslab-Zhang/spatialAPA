# What is spatialAPAdb
SpatialAPAdb ((http://www.biomedical-web.com/spatialAPAdb/home)) is a comprehensive atlas for exploring the spatial diversity and dynamics of alternative polyadenylation across across 18 human tissue organs and 76 disease states. SpatialAPAdb aims to explore the spatial diversity and dynamics of alternative polyadenylation (APA), polyA site index (PSI), weighted 3'-UTR length (WUL), and genes through analysis of 363 human spatial transcriptomics datasets from 56 projects, totaling 804,276 spots.

Figure 1. Workflow of spatialAPAdb.

# Source code of spatialAPAdb data processing.
We ingetrage gene expression matrix, APA expressiong matrix, PSI matrix and weighted 3'UTR length matrix for spatial transcriptome data (10x visium platform) using snakemake.
User can run this pipeline after all environment are pre-install for their raw data (fastq/fq with attachment slice).

## An integrated pipeline enable gene expression, APA/PSI, weighted 3'UTR length for 10x visium data
SCAPE (https://github.com/LuChenLab/SCAPE) should be pre-installed before using this pipeline.
Other dependancy of environment for running pipeline of spatialAPAdb are list in environment.yml file.

## Usage
```
git clone https://github.com/Omicslab-Zhang/spatialAPA.git
cd spatialAPA
conda env create -f environment.yml
snakemake -s /path/to /spatialAPA/snakefile/run_scape.pipeline.py \
  -j 80 -k --scheduler greedy --rerun-incomplete > /path/to/log/file/spatialAPA.log 2>&1
```

## Parameters
-	**-s /path/to/spatialAPA/snakefile/run_scape.pipeline.py:**  
Specifies the Snakemake workflow file (Snakefile) that contains the rules and steps for the pipeline. Replace /path/to/spatialAPA/snakefile/run_scape.pipeline.py with the actual path to your pipeline file.
-	**-j 80:**  
Defines the maximum number of jobs to run simultaneously. In this example, 80 threads are used.
-	**-k:**  
Keeps running independent jobs even if one fails, allowing partial progress in case of errors.
-	**--scheduler greedy:**  
Uses the greedy job scheduler to maximize the use of available resources during the workflow execution.
-	**--rerun-incomplete:**  
Ensures that any jobs that were interrupted or only partially completed in previous runs are rerun.
-	**/path/to/log/file/spatialAPA.log 2>&1:**  
Redirects both the standard output and error logs of the Snakemake execution to a log file for easy debugging and tracking. Replace /path/to/log/file/spatialAPA.log with the desired log file path.
Other options see [snakemake official documents](https://snakemake.readthedocs.io/en/v7.32.3/).

## Downstream analysis
The output of spatialAPA including gene expression/APA/PSI/WUL matrix are integrated  into a Seurat object in R, which can be loaded for downstream spatial transcriptomics analyses offered by Seurat.

## Further Information
To customize or debug the pipeline:
Modify config.yaml to set appropriate paths, parameters, and sample-specific details.
Check Snakemake's official documentation for advanced options and troubleshooting.
Ensure necessary dependencies and paths are set correctly, as defined in the run_scape.pipeline.py.
For more help, inspect specific rules or contact the pipeline maintainer.
