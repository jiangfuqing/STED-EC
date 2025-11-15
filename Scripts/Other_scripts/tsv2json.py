#!/usr/bin/env python
# coding: utf-8

# ###  generate metadata files needed for Seurat workflow for spatial datasets

import xml.etree.ElementTree as ET
import re
import pandas as pd
import argparse
import re

ap = argparse.ArgumentParser()
ap.add_argument("-i", "--sampleId", required=True, help="input_image")

args = vars(ap.parse_args())

sampleId = args["sampleId"]

#run in python39

#tsv to svg to json

#tsv to svg file
import os
import subprocess
import pickle

# use perl scripts
perl_script = "/public/home/fqjiang_gdl/SPARC/MISAR_to_Visium/Scripts/Other_scripts/tsv2svg.pl"

file_in = sampleId + '.tsv'
#output = subprocess.check_output(["perl", perl_script, "arg1", "arg2"])
output = subprocess.check_output(["/public/software/genomics/unstable/seeksoultools.1.2.0/bin/perl", perl_script, file_in ])

file_path = sampleId + '.svg'
with open(file_path, 'wb+') as file:
    file.write(output)


#svg to json

# Read image in svg format
dir_name = './'
svg_file = dir_name + sampleId + '.svg'
svg = ET.parse(svg_file)
root = svg.getroot()

# Prepare tissue_positions_list.csv file that compatible with Seurat
spots = pd.DataFrame(columns=['in_tissue', 'array_row', 'array_column', 
                              'pxl_col_in_fullres', 'pxl_row_in_fullres'])
col_reversed = False

spot_width = 0
spot_height = 0

for i, g in enumerate(root.findall("{http://www.w3.org/2000/svg}rect")):
    # print(g.attrib)
    x = float(g.attrib['x'])
    y = float(g.attrib['y'])
    spot_width = float(g.attrib['width'])
    spot_height = float(g.attrib['height'])
    
    if(g.attrib['fill'] == 'red'):
        in_tissue = 1
    else:
        continue
    #   in_tissue = 1
    
    # row = int(i/50)+1
    # col = int(i%50)+1
    col=int((x-100)/20+1)
    row=int((y-100)/20+1)
    if(col_reversed == True):
        col = 50-col+1
    
    x_c = x + spot_width/2
    y_c = y + spot_height/2
    
    spots = pd.concat([spots, pd.DataFrame.from_records([{'in_tissue': in_tissue, 'array_row': row-1, 'array_column': col-1, 'pxl_col_in_fullres': int(round(y_c)), 'pxl_row_in_fullres': int(round(x_c))}])])

spots.index = (spots['array_column']+1).astype(int).astype(str) + 'x' + (spots['array_row']+1).astype(int).astype(str)


# Read barcode.txt
barcode_txt = '/public/home/fqjiang_gdl/SPARC/barcode_file/SPARC_barcode_filter.csv'
barcodes = pd.read_csv(barcode_txt, sep=',')


# barcodes.columns = ['barcode', 'x', 'y']
# barcodes.index = barcodes['x'].astype(str) + 'x' + barcodes['y'].astype(str)
# barcodes = barcodes.drop(columns=['x', 'y'])


barcodes.index = barcodes['array'].astype(str)
barcodes = barcodes.drop(columns=['array'])


spots = pd.concat([barcodes, spots], axis=1)


spots = spots.dropna()


spots.to_csv(dir_name + '/tissue_positions_list.csv', index=False, header=False)


# Generate scalefactors_json.json file

import json

scalefactors = {"spot_diameter_fullres": spot_width, 
                "tissue_hires_scalef": 1.0, 
                "fiducial_diameter_fullres": spot_width, 
                "tissue_lowres_scalef": 1.0}

with open(dir_name + '/scalefactors_json.json', 'w') as outfile:
    json.dump(scalefactors, outfile)






