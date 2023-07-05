#!/usr/bin/env python
import argparse
import os
import pandas as pd

# import csv
# import numpy
# import shlex
# import shutil
# import subprocess
# import sys

parser = argparse.ArgumentParser(description=
    'Python script for parsing the bulk text file extracted from IPA. Example usage: python parseIPA.py /users/janedoe/sources/ipa_case.vs.control.txt /users/janedoe/output_folder')
parser.add_argument('--version', action='version', version='%(prog)s 0.1')
parser.add_argument('input_file', action= 'store', metavar='input_file', help='Enter the name and filepath of the bulk text file prefix name for data. Example: /users/janedoe/sources/ipa_case.vs.control.txt')
parser.add_argument('--outputfolder', action= 'store', dest='outputfolder', metavar='outputfolder', help='List folder to contain results. Example: /outputfolder')
parser.add_argument('--outputprefix', action= 'store', dest='outputprefix', metavar='outputprefix', help='Prefix for output files. By default, the prefix used for output is derived from the input file  Example: /outputfolder')
parser.add_argument('--basic', action= 'store', dest='basic_only', metavar='False', default= False, help='Default is FALSE, in which case only the tables for Upstream Regulators,Causal Networks, and Canonical Pathways will be saved. Specify TRUE if you wish to save the parsed tables from the entire input file')
parser.add_argument('--delim', action= 'store', dest='delim', metavar='delimiter', default= '\t', help='Use this delimiter for parsing the columns of results tables. It is highly recommended to not modify this default parameter.')
args = parser.parse_args()

# parameters
input_file = args.input_file
outputfolder = args.outputfolder
outputfile_prefix = args.outputprefix
basic_only = args.basic_only
delim = args.delim


# input_file ="C:\\Users\\Molly\\Documents\\erlich\\Microglia_May2023\\ipa\\ipa_mor.vs.sal_pnd14.txt"
## remove the file path to obtain the file prefix
if outputfile_prefix is None: 
    outputfile_prefix=input_file.split(os.path.sep)[-1]
else:
    print("Using the file_prefix: "+outputfile_prefix)

## Output directory is generated from the input filepath if not explictly provided
if outputfolder is None: 
    outputfolder=os.path.sep.join(input_file.split(os.path.sep)[:-1])
else:
    print("Using the output directory: "+outputfolder)
## Derive the file prefix to be used later in saving the files
out_prefix=outputfolder + os.path.sep + outputfile_prefix

with open(input_file) as f:
    lines=f.read()

curdict=dict()
b=lines.split("\n\n")
for i in range(len(b)):
    # print(i)
    a=b[i].split("\n")
    subdict=dict()
    for i2,line in enumerate(a):
        a2=line.split("\t")
        subdict[i2]=a2
    curdict[i]=pd.DataFrame.from_dict(subdict,orient='index')

re_dict=dict()
for i in range(len(b)):
    # print(i)
    if "Upstream Regulators for My Projects" in curdict[i].loc[0,0] :
        re_dict["upreg"]=curdict[i].loc[1:,:]
    if "Causal Networks for My Projects" in curdict[i].loc[0,0] :
        re_dict["causal"]=curdict[i].loc[1:,:]
    if "Canonical Pathways for My Projects" in curdict[i].loc[0,0] :
        re_dict["canonical"]=curdict[i].loc[1:,:]
    if "Diseases and Bio Functions for My Projects" in curdict[i].loc[0,0] :
        re_dict["biofun"]=curdict[i].loc[1:,:]
    if "Tox Functions for My Projects" in curdict[i].loc[0,0] :
        re_dict["toxfun"]=curdict[i].loc[1:,:]
    if "Regulator Effects for My Projects" in curdict[i].loc[0,0] :
        re_dict["regeff"]=curdict[i].loc[1:,:]
    if "Analysis Ready Molecules for My Projects" in curdict[i].loc[0,0] :
        re_dict["amols"]=curdict[i].loc[1:,:]

## We need to clean up the column names, reassign as the first column
for key,value in re_dict.items():
    print(key)
    print(value.empty)
    if not value.empty:      
        value.columns=value.iloc[0,:]
        value=value.iloc[1:,:]
        re_dict[key]=value

re_dict["upreg"].to_csv(out_prefix+"_upstreamregulators.txt",sep="\t",index=False)
re_dict["causal"].to_csv(out_prefix+"_causalnetworks.txt",sep="\t",index=False)
re_dict["canonical"].to_csv(out_prefix+"_canonicalpathways.txt",sep="\t",index=False)
if not basic_only:
    re_dict["biofun"].to_csv(out_prefix+"_diseases.txt",sep="\t",index=False)
    re_dict["toxfun"].to_csv(out_prefix+"_toxfunctions.txt",sep="\t",index=False)
    re_dict["regeff"].to_csv(out_prefix+"_regulatoreffects.txt",sep="\t",index=False)
    re_dict["amols"].to_csv(out_prefix+"_analysisready.txt",sep="\t",index=False)

## Let the user know if they are generating all files or simply the basic set
# if basic_only:
#     print("Only basic files will be saved.")
# else:
#     print("All files will be saved.")
