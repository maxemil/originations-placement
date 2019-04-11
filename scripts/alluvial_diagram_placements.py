#!/usr/bin/env python3
"""
  19-04-04, Max Schön, <max-emil.schon@icm.uu.se>
  example usage:
"""
import os
import argparse
import pandas as pd
import rpy2.robjects as robjects
import rpy2.robjects.packages as rpackages

parser = argparse.ArgumentParser()

parser.add_argument("placements", type=str,nargs='+',
                    help="gappa csv output files")
parser.add_argument("-a", "--annotations" type=str,
                    help="EggNOG OG annotation file")

args = parser.parse_args()


annot = {line.split('\t')[0]:eval(line.split('\t')[4]) for line in open(args.annotations)}

tax = []
for placement in args.placements:
    base = os.path.basename(placement).split('.')[0]
    df = pd.read_csv(placement, sep='\t')
    tax.append(df['taxopath'][df['fract'].argmax()].split(';')[-1])


tit <- as.data.frame(Titanic, stringsAsFactors = FALSE)
alluvial(tit[,1:4], freq=tit$Freq,
         col = ifelse(tit$Survived == "Yes", "orange", "grey"),
         border = ifelse(tit$Survived == "Yes", "orange", "grey"),
         hide = tit$Freq == 0,
         cex = 0.7
)
