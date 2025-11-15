import diopy
import argparse
import re

ap = argparse.ArgumentParser()
ap.add_argument("-i", "--sampleId", required=True, help="input_image")

args = vars(ap.parse_args())

sampleId = args["sampleId"]

# read data
adata_omics1 = diopy.input.read_h5(file = sampleId + '-RNA.h5', assay_name='spatial')

adata_omics2 = diopy.input.read_h5(file = sampleId + '-ATAC.h5', assay_name='spatial')

adata_omics3 = diopy.input.read_h5(file = sampleId + '-GeneScore.h5', assay_name='spatial')

# save data
adata_omics1.write(sampleId + '-RNA.h5ad')

adata_omics2.write(sampleId + '-ATAC.h5ad')

adata_omics3.write(sampleId + '-GeneScore.h5ad')