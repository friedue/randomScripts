#!/bin/sh
##############################
# (c) FrieDue, Aug 2015
##############################
# This script turns the wig file of eRRBS runs produced by the Epicore into bigWig files.
# To achieve this, the wig file is first split into its 4 parts:
# methCpG, methCpG_R, NOmethCpG, NOmethCpG_R - these are based on the following
# track lines found within the EpicCores Wig:
#	track type=wiggle_0 name="B1 methCpG" visibility=full color=190,0,0 autoScale=off viewLimits=0.0:100.0 yLineMark=0.0 yLineOnOff=on maxHeightPixels=75:75:20 priority=01
#	track type=wiggle_0 name="B1 methCpG_R" visibility=full color=150,20,20 autoScale=off viewLimits=0.0:-100.0 yLineMark=0.0 yLineOnOff=on maxHeightPixels=75:75:20 priority=01
#	track type=wiggle_0 name="B1 NOmethCpG" visibility=full color=0,0,190 autoScale=off viewLimits=0.0:100.0  yLineMark=0.0 yLineOnOff=on maxHeightPixels=75:75:20 priority=02
#	track type=wiggle_0 name="B1 NOmethCpG_R" visibility=full color=20,20,150 autoScale=off viewLimits=0.0:-100.0  yLineMark=0.0 yLineOnOff=on maxHeightPixels=75:75:20 priority=02
# MethCpG and methCpG_R values are combined into one single wig file, which is subsequently
# converted into a bigWig using UCSC tools.
##############################
# usage: EpiWig_to_UCSCbigWig.sh B1_CpG.wig hg19.chromInfo.txt convertMinus
##########################
EpicoreWig=${1} # Epicore eRRBS wig file
ChromInfo=${2} # for bigWig conversion
##########################

# split Epicore wig file into individual files per track
# using the information from name="xxxx" for the naming of the output files (stored in $2)
echo "Splitting ${EpicoreWig} into 4 tracks"
grep -v browser ${EpicoreWig} | awk 'BEGIN{FS="\""} /track/ {gsub(" ","_",$2);out="tmp."$2}{OFS="\"";print > out}'

rm tmp*NOmeth*

if [ "$3" = convertMinus ]; then
	TMP=`ls tmp.*methCpG_R | sed 's/.*tmp\.\(.*\_methCpG\_R\)/\1/g'`
	echo "converting negative strand CpG-me values into positive integers"
	awk '{for (i=1; i<=NF; i++) if ($i <= 0) $i = -$i; print}' tmp.*methCpG_R > tmp2.${TMP}
	rm tmp.*methCpG_R
	mv tmp2.${TMP} tmp.${TMP}
fi

# to merge fwd and rev strand, I first split the files per chromosome
# (indicated after "variableStep") and merge them chromosome by chromosome
for i in methCpG methCpG_R
do
File=tmp*${i}
Sample=`echo ${File} | sed 's/tmp\.\([a-zA-Z0-9_]*\)/\1/'`
echo "Splitting chromosomes for ${Sample}"
grep -v track ${File} | awk -v var1=${Sample} 'BEGIN {FS="="}/variableStep/{x="tmp."var1"_"$2}{print > x}'
done

# merge fwd and rev
Sample_short=`echo ${Sample} | sed 's/_methCpG.*//g'`
echo "merging forward and reverse strand counts for ${Sample_short}"
for CHR in `ls tmp.*methCpG_chr* | awk 'BEGIN {FS="_"}{print $NF}'`
do
cat tmp.${Sample_short}_methCpG_${CHR} tmp.${Sample_short}_methCpG_R_${CHR} | sort -k1,1n | uniq >> tmp.${Sample_short}_freqC.wig
done

# turn wig into bigwig
echo "generating bigWig of cytosine counts: ${Sample_short}_freqC.bw"
wigToBigWig tmp.${Sample_short}_freqC.wig ${ChromInfo} ${Sample_short}_freqC.bw

# cleaning up
rm tmp.${Sample_short}*

