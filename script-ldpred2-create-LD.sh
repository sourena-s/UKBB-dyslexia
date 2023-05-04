#!/bin/bash

module unload R
module load R
Rscript script-ldpred2-create-LD.R $1
