#!/bin/bash

# 输入参数
GSM=$1
FASTQS=$2
TRANSCRIPTOME=$3
CORES=$4
MEM=$5
SLIDE=$6
OUTDIR=$7
IMAGE=$8
DARKIMAGE=$9
LOUPE=$10

# 运行 spaceranger count
spaceranger count \
    --id=${GSM} \
    --fastqs=${FASTQS} \
    --sample=${GSM} \
    --transcriptome=${TRANSCRIPTOME} \
    --localcores=${CORES} \
    --localmem=${MEM} \
    ${SLIDE} \
    ${OUTDIR} \
    ${IMAGE:+--image ${IMAGE}} \
    ${DARKIMAGE:+--darkimage ${DARKIMAGE}} \
    ${LOUPE:+--loupe-alignment ${LOUPE}}
