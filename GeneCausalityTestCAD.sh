#!/bin/bash

WORKDIR=~/CADcausalitytest
GENENAME="$(cat $1)"
CHR=$2
DISTANCE=$3
TRESHOLD=$4
NELSON=UKBB.GWAS1KG.EXOME.CAD.SOFT.META.PublicRelease.300517.txt
CC4D=~/CADcausalitytest/CARDIOGRAMC4D


if [ ! -d $WORKDIR ]
then
mkdir $WORKDIR
fi

if [ ! -d $CC4D ]
then
mkdir $CC4D
fi

cd $CC4D

if [ ! -f $NELSON ]
then
wget 
echo "Unpacking CARDIOGRAM plus C4D data, Nelson et al."
unzip ...
fi


awk '{if ($3=="'"$CHR"'") print $0}' UKBB.GWAS1KG.EXOME.CAD.SOFT.META.PublicRelease.300517.txt > UKBB.GWAS1KG.EXOME.CAD.SOFT.META.PublicRelease.300517.chr$CHR.txt

#write R script to get ENSEMBL id, needs biomaRt in R

echo "#!/usr/bin/Rscript
library(biomaRt)
listMarts(host=\"grch37.ensembl.org\")
ensembl = useMart(\"ENSEMBL_MART_ENSEMBL\",dataset=\"hsapiens_gene_ensembl\", host=\"grch37.ensembl.org\")
id_merge = getBM(attributes=c(\"ensembl_gene_id\",\"external_gene_name\",\"transcript_start\",\"transcript_end\"),mart=ensembl)
write.table(id_merge, file=\"id_merge.txt\", sep = \"\t\", quote =F, col.names=F, row.names=F)
" > script.r

#run R script

chmod 775 script.r
./script.r

grep $GENENAME id_merge.txt > gene_id_merge

START="$(awk 'BEGIN{min=100000000000000000}{if($3<min){min=$3;out=$3}}END{print out}' gene_id_merge)"
END="$(awk 'BEGIN{max=0}{if($4>max){max=$4;out=$4}}END{print out}' gene_id_merge)"

echo "Composite gene start:" $START
echo "Composite gene end:" $END
echo "Extending with the distance:" $DISTANCE "bp"

awk 'NR==FNR {if ($4>="'"$START"'"-"'"$DISTANCE"'" && $4<="'"$END"'"+"'"$DISTANCE"'") print $0}' UKBB.GWAS1KG.EXOME.CAD.SOFT.META.PublicRelease.300517.chr$CHR.txt  > chr$CHR.$GENENAME.region.txt
 
cut -f2,5,10 chr$CHR.$GENENAME.region.txt > SNP_effect.alele_pval.txt

echo "Selelcting p-value threshold" $TRESHOLD
awk -v var="$TRESHOLD" '{if ($3<=var) print $0}' SNP_effect.alele_pval.txt > SNP_effect.alele_pval.threshold.txt

#create genotype files

cd $WORKDIR/HCASMC_genotypes/vcf

        grep CHROM phased_and_imputed.chr$CHR.vcf > HEADER.txt
        awk '{$1=$2=$3=$4=$5=$6=$7=$8=$9=""; print $0}' HEADER.txt > HEADER.txt.cut
	#cat HEADER.txt.cut
while read -r a b c; do
        grep "$a	" phased_and_imputed.chr$CHR.vcf > SNP.txt
	#cat SNP.txt
        #clean genotype files
        sed -i -E 's/:[0-9.]*:[0-9.]*,[0-9.]*,[0-9.]*//g' SNP.txt 
	#cat SNP.txt
       
#greb reference and alternative aleles
        REF="$(awk '{printf $4}' SNP.txt)"
        ALT="$(awk '{printf $5}' SNP.txt)"
	#echo $REF
	#echo $ALT

        EFFECT=$b
        #echo $EFFECT
	
	sed -i "s/0|0/$REF$REF/g" SNP.txt
        sed -i "s/0|1/$REF$ALT/g" SNP.txt
        sed -i "s/1|0/$REF$ALT/g" SNP.txt
        sed -i "s/1|1/$ALT$ALT/g" SNP.txt
	#cat SNP.txt
        
        awk '{$1=$2=$3=$4=$5=$6=$7=$8=$9=""; print $0}' SNP.txt > SNP.txt.cut 
        #cat SNP.txt.cut

	sed -i "s/"$EFFECT"/1/g" SNP.txt.cut
        sed -i -E "s/[ATGC]/0/g" SNP.txt.cut
        #cat SNP.txt.cut

        cat HEADER.txt.cut SNP.txt.cut > GENOTYPES.$a.txt
	#cat GENOTYPES.$a.txt

done < $CC4D/SNP_effect.alele_pval.threshold.txt

cat GENOTYPES.rs* > GENOTYPES.combined
awk 'NR%2==0' GENOTYPES.combined > GENOTYPES.combined.even
cat HEADER.txt.cut GENOTYPES.combined.even > GENOTYPES.combined.even.HEADER

sed -i 's/11/2/g' GENOTYPES.combined.even.HEADER
sed -i -E 's/(10|01)/1/g' GENOTYPES.combined.even.HEADER
sed -i 's/00/0/g' GENOTYPES.combined.even.HEADER

####awk 'BEGIN{print "count", "lineNum"}{print gsub(/t/,"") "\t" NR}' file
####awk -F'|' -v fld=2 'BEGIN{print "count", "lineNum"}{print gsub(/t/,"",$fld) "\t" NR}' file

echo "#!/usr/bin/Rscript
library(\"ggplot2\")
data<-read.table (file=\"GENOTYPES.combined.even.HEADER\", head=T, check.names=F)
data.tr<-t(data)
data.sum<-rowSums(data.tr)







