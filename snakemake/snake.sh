#!/usr/bin/bash
nohup snakemake --cluster "qsub -q all" -j 64 --snakefile snake.py >snake.log &
