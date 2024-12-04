#this is a script to get the top hit of the blast result blasted againt the database of interest
#_______________________________________________________________________________________________#
import os
import pandas as pd
import argparse
#_______________________________________________________________________________________________#
#data_file
#blst_hits = "/mnt/gs21/scratch/ranawee1/Potato_colab/genomes/finding_homology/VER_blast_against_DM1S1.txt"

#adding parser to the script
parser = argparse.ArgumentParser(description='Get the top hit of the blast result')
parser.add_argument('-i', '--input', help='blast result file', required=True)
args = parser.parse_args()
blst_hits = args.input

#open the blast result file
df_blst_hits = pd.read_csv(blst_hits, sep='\t', header=None)

#get the top hit for unique egenes in column 0
df_top_hit = df_blst_hits.groupby(0).first()
#save the first and second columns and 11th column to a txt file
df_top_hit[[1,10]].to_csv("DM1S1_vs_VER_blast_top_hit.txt", sep='\t', header=None)



