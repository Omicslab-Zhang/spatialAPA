# What is the spatial APA data analysis pipeline?
To quantify and integrate spatial APA data with gene expression, we developed a comprehensive spatial APA data analysis pipeline, called spatialAPA, built upon Space Ranger, spatialAPA, and Snakemake. This pipeline encompasses multiple steps, including:
(1) Quality control of raw sequencing data (FASTQ/FQ files with corresponding tissue slice images).
(2) Quantification of APA events and gene expression.
(3) Calculation of Polyadenylation Site Index (PSI) and Weighted 3' UTR Length (WUL) scores for APA events.
(4) Downstream analysis: including cluster, spatial variable feature analysis, and GO/KEGG enrichment analysis.

Users can easily run spatialAPA after configuring the required environment for their raw data ( fastq/fq of 10x visium spatial transcripomics data with attachment slice). The outputs, including APA matrices, PSI, WUL, and gene expression data, are seamlessly integrated into a Seurat object in R. This integration enables downstream spatial transcriptomics analyses using Seurat's extensive toolkit.

Using the spatialAPA pipeline, we constructed a comprehensive spatially resolved atlas of APA, analyzing 804,276 transcriptomic spots from 363 sections, spanning 56 projects, 18 organs, and 76 disease states. The resulting atlas provides unprecedented insights into the spatial dynamics of APA and is available for exploration at [SpatialAPAdb](http://www.biomedical-web.com/spatialAPAdb/home).

![image](https://github.com/Omicslab-Zhang/spatialAPA/blob/main/STAPADB_flow.drawio.png)  
<p align="center">
  <strong>Figure 1. Workflow of spatialAPAdb.</strong>
</p>

## Usage of the spatialAPA pipeline
```
git clone https://github.com/Omicslab-Zhang/spatialAPA.git
cd spatialAPA
conda env create -f environment.yml
snakemake -s /path/to/spatialAPA/snakefile/run_spatialAPA.pipeline.py \
  -j 80 -k --scheduler greedy --rerun-incomplete > /path/to/log/file/run_spatialAPA.log 2>&1
```

## Explanation of Key Command-Line Options for Running the Spatial APA Pipeline
- **-s /path/to/spatialAPA/snakefile/run_spatialAPA.pipeline.py:**  
Specifies the Snakemake workflow file (Snakefile) containing the rules and steps for the pipeline. Replace /path/to/spatialAPA/snakefile/run_spatialAPA.pipeline.py with the actual file path of your pipeline.
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
The [SCAPE](https://github.com/LuChenLab/SCAPE) tool should be pre-installed before using this pipeline. Other dependancy of environment for running pipeline of spatialAPAdb are list in environment.yml file.
For additional options and detailed instructions, refer to the [Snakemake documentation](https://snakemake.readthedocs.io/en/v7.32.4/) or the pipeline's user guide.

## Further Information
To customize or debug the pipeline:
- Modify config.yaml to set appropriate paths, parameters, and sample-specific details.
- Check Snakemake's official documentation for advanced options and troubleshooting.
- Ensure necessary dependencies and paths are set correctly, as defined in the run_spatialAPA.pipeline.py.
- For more help, inspect specific rules or contact the pipeline maintainer.
