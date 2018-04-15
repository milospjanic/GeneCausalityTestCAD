# GeneCausalityTestCAD

GeneCausalityTest for CAD (coronary artery disease) is a combined bash/awk/R script for defining causality of a gene for a given trait, in this case CAD, given the directionality of change of expression level with the increasing number of risk GWAS SNPs. Script uses CAD GWAS data from the most recent meta analysis (Nelson et al. Association analyses based on false discovery rate implicate new loci for coronary artery disease. Nat Genet 2017. doi: 10.1038/ng.3913) to select risk SNPs and define risk alleles and HCASMC eQTL data to perform regression analysis.

User provides as arguments:

1. gene of interest
2. chromosome where the gene is located
3. defines p-value threshold to select SNPs from Nelson et al. meta analysis for CAD 
4. defines plus/minus region to select SNPs from, starting from composite gene coordinates obtained after collapsing all gene isoforms into a composite gene.

# Usage

To run the script download the .sh file
<pre>
wget https://raw.githubusercontent.com/milospjanic/GeneCausalityTestCAD/master/GeneCausalityTestCAD.sh
chmod 755 HCASMCeQTLviewer.sh
</pre>

Place the script in your home or any other folder. The script will create ~/CADcausalitytest as its working folder, and three subfolders: CARDIOGRAMC4D, HCASMC_expr and HCASMC_genotypes. HCASMC_expr will contain per-gene RNAseq read counts for each HCASMC sample, while HCASMC_genotypes/vcf/ contains whole genome sequencing vcf files of HCASMC samples. CARDIOGRAMC4D folder will contain summary data from Nelson et al. Script will check if all three folders and data sets are present and if not download.

# Examples



 
