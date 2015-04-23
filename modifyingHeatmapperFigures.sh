# crops the figures that are output from heatmapper (default setting, "heatmapper only") and adds a dashed line in the center of the remaining figure
# also adds a title that should be specified
PICTURE=${1}
TITLE=${2}
declare OUTNAME=${PICTURE%\.png*}
convert ${PICTURE} -crop 265x1030+22+25 -trim +repage -fill none -stroke black -strokewidth 1.5 -draw "stroke-dasharray 5 5 line 105,982 105,0" ${OUTNAME}.trimmed.png
montage -title ${TITLE} ${OUTNAME}.trimmed.png -geometry +0.1+0.1 ${OUTNAME}.modified.png
rm ${OUTNAME}.trimmed.png
