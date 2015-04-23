#!/bin/sh
PEAKS=ClusterD.bed
APPENDIX=50bp_ClusterD
SIGNALTYPE=log2ratio
SIGNALMATRIX=${SIGNALTYPE}_${APPENDIX}
DISTANCE=euclidean
CLUSTER=ward

#1. signalMatrix generation mlES
echo "producing signal matrix for ${CELL}"
for CHIP in MOF_ES H3K27ac_ES H3K4me2_ES H3K4me1_ES
do
SIGNAL=${CHIP}.bw
/package/deepTools-1.5/bin/computeMatrix scale-regions -S ${SIGNAL} -R ${PEAKS} -out compMatrix_${SIGNAL}.gz --outFileNameMatrix matrix_${CHIP}_${SIGNALMATRIX}.tab -bs 50 -m 1400 -a 0 -b 0 --missingDataAsZero
echo matrix_${CHIP}_${SIGNALMATRIX}.tab >> filelist
done

#2. signal Matrix rank transformation
echo "performing signal matrix rank transformation"
/data/projects/muehlpfordt/scripts/Cluster0_matrixProcessing.R filelist ${SIGNALTYPE}_${APPENDIX}

#3. hclust
echo "performing hclust"
/data/projects/muehlpfordt/scripts/Cluster1_computeHierarchicalClustering_FM_2013-03-11.R ${SIGNALMATRIX}.tab ${DISTANCE} ${CLUSTER} ${SIGNALTYPE} new

#4. pruning
echo "performing dendrogram pruning"
for i in 3 4 5 6 7 8 9
do
/data/projects/muehlpfordt/scripts/Cluster2a_split_in_groups_plot_dendrogram_FM_outputUnsorted.R ${SIGNALMATRIX}.tab hclust_${SIGNALTYPE}_${DISTANCE}_${CLUSTER}.Robj ${PEAKS} ${i} clusters_${i}.bed dendrogram_${i}.pdf summary
/package/BEDTools/bin/bedtools shuffle -i ${PEAKS} -g /data/projects/misc/genomes/Mm/mm9/sizes/mm9_no_random.genome | head -2500 >> clusters_${i}.bed 
for CHIP in MOF_ES H3K27ac_ES H3K4me2_ES H3K4me1_ES
do
SIGNAL=${CHIP}.bw
/package/deepTools-1.5/bin/computeMatrix reference-point -bs 50 -a 2500 -b 2500 -S ${SIGNAL} -R clusters_${i}.bed --referencePoint center --missingDataAsZero -out compMatrix_${i}Clusters_${CHIP}.gz
/package/deepTools-1.5/bin/heatmapper -m compMatrix_${i}Clusters_${CHIP}.gz -T "${i} Clusters: ${CHIP} signal" -out ${i}Clusters_${CHIP}.png --sortRegions descend --whatToShow "heatmap and colorbar"
done
montage ${i}Clusters_MOF_ES.png ${i}Clusters_H3K27ac_ES.png ${i}Clusters_H3K4me2_ES.png ${i}Clusters_H3K4me1_ES.png  -tile x1 -geometry +0.1+0.1 ${i}Clusters_MOFetc_${SIGNALTYPE}.png
done

for i in 3 4 5 6 7 8 9
do
rm ${i}Clusters_MOF_ES.png ${i}Clusters_H3K27ac_ES.png ${i}Clusters_H3K4me2_ES.png ${i}Clusters_H3K4me1_ES.png
done

montage 3Clusters_MOFetc_${SIGNALTYPE}.png 4Clusters_MOFetc_${SIGNALTYPE}.png 5Clusters_MOFetc_${SIGNALTYPE}.png  6Clusters_MOFetc_${SIGNALTYPE}.png 7Clusters_MOFetc_${SIGNALTYPE}.png 8Clusters_MOFetc_${SIGNALTYPE}.png 9Clusters_MOFetc_${SIGNALTYPE}.png -tile x7 -geometry +0.1+0.1 AllClusters_MOFetc_${SIGNALTYPE}.png
