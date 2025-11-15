#!/usr/bin/bash
#MISAR_to_Visium
#Fuqing Jiang

sampleId=$1
Species=$2

cd Data

#counts to seurat object
python ~/Scripts/Figure_filter.py -i $sampleId
Rscript ~/Scripts/counts2tsv.R $sampleId
python ~/Scripts/tsv2json.py -i $sampleId
python ~/Scripts/resize_img.py -in ${sampleId}-HE.jpg -out ${sampleId}-HE.png
rm $sampleId.tsv $sampleId.svg
#generate sampleId-RNA.h5 sampleId-ATAC.h5 sampleId-GeneScore.h5
Rscript ~/Scripts/Visium_Seurat_GeneScore.R $sampleId $Species
#transform h5 to h5ad file using diorpy in python3.7
#since the numpy version needed by diorpy was 1.20.0, but SpatialGlue need numpy version was 1.24.0
python ~/Scripts/diopy_h5.py -i $sampleId
rm $sampleId-RNA.h5 $sampleId-ATAC.h5 $sampleId-GeneScore.h5

#run SpatialGlue
python $Main/Scripts/Other_scripts/SpatialGlue_MISAR.py -i $sampleId
