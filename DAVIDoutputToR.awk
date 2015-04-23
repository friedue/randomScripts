BEGIN{
    CLUST=0
    }
/^$/ {next} # makes awk ignore blank lines
{
    if($1=="Annotation" && $3==(CLUST+1)){
         CLUST=$3; ENRICH=$6; next
        }
    else
        if(/Category/ ) {next} # skips lines with Category or Enrichment
        else
	OFS="\t"
        print "Cluster"CLUST,ENRICH,$0
}
