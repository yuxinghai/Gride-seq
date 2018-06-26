#!/usr/bin/bash
nohup snakemake -j 64 --snakefile snake.py >snake.log &
