#!/bin/sh 

###20190321 -- CAFE (Correlated Allele Frequency Evolution) (?)

#20190321
#TOC

 - 20190321: Vs1


#20181106
#Re-built website structure/platform using `workflowr` (see: https://jdblischak.github.io/workflowr/index.html)
R -q -e "library(\"workflowr\"); wflow_start(\"/users/mturchin/LabMisc/RamachandranLab/InterPath\", existing=TRUE, git=FALSE); wflow_build();"
#Useful workflowr related commans
R -q -e "library(\"workflowr\"); wflow_build();"
#scp -p /users/mturchin/LabMisc/RamachandranLab/InterPath/docs/* /Users/mturchin20/Documents/Work/LabMisc/RamachandranLab/InterPath/docs/.
#From MacBook Pro
#rsync -zvhia mturchin@ssh.ccv.brown.edu:/users/mturchin/LabMisc/RamachandranLab/InterPath/docs/. /Users/mturchin20/Documents/Work/LabMisc/RamachandranLab/InterPath/docs/.  

#20190321
#Conda environment -- using `MultiEthnicGWAS`

#20200114
#Updated naming convention from `CAFE` to `ShiftAE`


#20190321
#Vs1

#Idea -- SNPs undergoing drift should randomly fluctuate to both higher and lower frequencies 'equally' (modulo population demographic misc) from their ancestral state over time. SNPs undergoing weak selective pressures in a given direction should experience less overall variance in their change in frequency since their ancestral state (alleles that are favored get pushed up; if you don't know orientation properly, can just use absolute difference since amount pushed up or amount pushed down are equivalent; looking at overall standard difference amongst 'AF change since ancestral state'; AF affects how much has changed since ancestral, but will be doing matched permutations on ancestral AF so within the final score, that's fine; the idea is that all alleles have similar changes in AF since ancestral state when they're being weakly selected, so they are more similar than random chance alone given a normal distribution variance)

mkdir /users/mturchin/LabMisc/RamachandranLab/InterPath/ShiftAE
mkdir /users/mturchin/LabMisc/RamachandranLab/InterPath/ShiftAE/Vs1
mkdir /users/mturchin/LabMisc/RamachandranLab/InterPath/ShiftAE/Vs1/Analyses

for pop1 in `echo "CEU GBR TSI IBS YRI ESN"`; do 
	echo $pop1

	if [ ! -d /users/mturchin/data/1000G/subsets/$pop1/AFs ] ; then
		mkdir /users/mturchin/data/1000G/subsets/$pop1/AFs
	fi
	
	for chr1 in `echo {1..22} X`; do

#		plink --bfile /users/mturchin/data/1000G/subsets/$pop1/${pop1}.chr${chr1}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.SNPs --freq --out /users/mturchin/data/1000G/subsets/$pop1/AFs/${pop1}.chr${chr1}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.SNPs  
		gzip -f /users/mturchin/data/1000G/subsets/$pop1/AFs/${pop1}.chr${chr1}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.SNPs.frq

	done
done

mkdir $HOME/data/mturchin/Neale2017/Vs2
mv $HOME/data/mturchin/Neale2017/* $HOME/data/mturchin/Neale2017/Vs2/.
#20200114 NOTE: below move to $HOME/data/mturchin/Data directory setup
#mkdir $HOME/data/mturchin/Data
#mv $HOME/data/mturchin/Neale2017 $HOME/data/mturchin/Data/.
#mv $HOME/data/mturchin/Loh2017 $HOME/data/mturchin/Data/.
#mkdir $HOME/data/mturchin/Data/GIANT2014_5
#mv /users/mturchin/Data2/GIANT/*vs1.lt* /users/mturchin/data/mturchin/Data/GIANT2014_5/.
#mv /users/mturchin/Data2/GIANT/GIANT2014_5.Height.edits.txt.gz /users/mturchin/data/mturchin/Data/GIANT2014_5/.
#cp -p /users/mturchin/Data2/GIANT/* /users/mturchin/data/mturchin/Data/GIANT2014_5/.

#From midway
#scp -p /project/mstephens/mturchin20/gluster/data/external_public_supp/GIANT2014_5/GIANT_HEIGHT_Wood_et_al_2014_publicrelease_HapMapCeuFreq.wUCSCGenomeBrowser_dbSNP130.vs1.txt.gz mturchin@ssh.ccv.brown.edu:/users/mturchin/Data2/GIANT/.

zcat $HOME/data/mturchin/Data/Neale2017/Vs2/50.assoc.tsv.gz | awk '{ if ($9 < 5e-9) { print $0 } } ' | gzip > $HOME/data/mturchin/Data/Neale2017/Vs2/50.assoc.lt5eNeg9.tsv.gz
zcat $HOME/data/mturchin/Data/Neale2017/Vs2/50.assoc.tsv.gz | awk '{ if ($9 < 1e-4) { print $0 } } ' | gzip > $HOME/data/mturchin/Data/Neale2017/Vs2/50.assoc.lt1eNeg4.tsv.gz
zcat $HOME/data/mturchin/Data/Loh2017/body_HEIGHTz.sumstats.gz | awk '{ if ($10 < 5e-9) { print $0 } } ' | gzip > $HOME/data/mturchin/Data/Loh2017/body_HEIGHTz.sumstats.lt5eNeg9.gz
zcat $HOME/data/mturchin/Data/Loh2017/body_HEIGHTz.sumstats.gz | awk '{ if ($10 < 1e-4) { print $0 } } ' | gzip > $HOME/data/mturchin/Data/Loh2017/body_HEIGHTz.sumstats.lt1eNeg4.gz
zcat $HOME/data/mturchin/Data/GIANT2014_5/GIANT_HEIGHT_Wood_et_al_2014_publicrelease_HapMapCeuFreq.wUCSCGenomeBrowser_dbSNP130.vs1.txt.gz | awk '{ if ($7 < 5e-9) { print $0 } } ' | gzip > $HOME/data/mturchin/Data/GIANT2014_5/GIANT_HEIGHT_Wood_et_al_2014_publicrelease_HapMapCeuFreq.wUCSCGenomeBrowser_dbSNP130.vs1.lt5eNeg9.txt.gz
zcat $HOME/data/mturchin/Data/GIANT2014_5/GIANT_HEIGHT_Wood_et_al_2014_publicrelease_HapMapCeuFreq.wUCSCGenomeBrowser_dbSNP130.vs1.txt.gz | awk '{ if ($7 < 5e-8) { print $0 } } ' | gzip > $HOME/data/mturchin/Data/GIANT2014_5/GIANT_HEIGHT_Wood_et_al_2014_publicrelease_HapMapCeuFreq.wUCSCGenomeBrowser_dbSNP130.vs1.lt5eNeg8.txt.gz
zcat $HOME/data/mturchin/Data/GIANT2014_5/GIANT_HEIGHT_Wood_et_al_2014_publicrelease_HapMapCeuFreq.wUCSCGenomeBrowser_dbSNP130.vs1.txt.gz | awk '{ if ($7 < 1e-4) { print $0 } } ' | gzip > $HOME/data/mturchin/Data/GIANT2014_5/GIANT_HEIGHT_Wood_et_al_2014_publicrelease_HapMapCeuFreq.wUCSCGenomeBrowser_dbSNP130.vs1.lt1eNeg4.txt.gz

##zcat $HOME/data/mturchin/Data/Neale2017/Vs2/50.assoc.tsv.gz | sed 's/:/ /g' | awk '{ print $5 "\t" $3 "\t" $4 "\t" $9 "\t" $10 "\t" $12 }' | grep -v pval | cat <(echo -e "rsid\ta1\ta2\tbeta\tse\tpval") - | gzip > $HOME/data/mturchin/Data/Neale2017/Vs2/50.assoc.edits1.tsv.gz
zcat $HOME/data/mturchin/Data/Neale2017/Vs2/50.assoc.tsv.gz | sed 's/:/ /g' | awk '{ print $5 "\t" $3 "\t" $4 "\t" $9 "\t" $10 "\t" $12 }' | grep -v pval | cat <(echo -e "rsid\ta1\ta2\tbeta\tse\tpval") - | gzip > $HOME/data/mturchin/Data/Neale2017/Vs2/Neale2017.Height.edits.txt.gz
ln -s $HOME/data/mturchin/Data/Neale2017/Vs2/Neale2017.Height.edits.txt.gz $HOME/data/mturchin/Data/Neale2017/Neale2017.Height.edits.txt.gz
zcat $HOME/data/mturchin/Data/Loh2017/body_HEIGHTz.sumstats.gz | awk '{ if ($4 != $6) { $8 = -1 * $8; val1 = $4; val2 = $5; $4 = val2; $5 = val1; }; print $1 "\t" $4 "\t" $5 "\t" $8 "\t" $9 "\t" $10 }' | grep -v SNP | cat <(echo -e "rsid\ta1\ta2\tbeta\tse\tpval") - | gzip > $HOME/data/mturchin/Data/Loh2017/Loh2017.Height.edits.txt.gz
zcat $HOME/data/mturchin/Data/GIANT2014_5/GIANT_HEIGHT_Wood_et_al_2014_publicrelease_HapMapCeuFreq.wUCSCGenomeBrowser_dbSNP130.vs1.txt.gz | awk '{ print $1 "\t" $2 "\t" $3 "\t" $5 "\t" $6 "\t" $7 }' | grep -v Allele | cat <(echo -e "rsid\ta1\ta2\tbeta\tse\tpval") - | gzip > $HOME/data/mturchin/Data/GIANT2014_5/GIANT2014_5.Height.edits.txt.gz

##zcat $HOME/data/mturchin/Data/Loh2017/body_HEIGHTz.sumstats.lt5eNeg9.gz | R -q -e "library(\"devtools\"); devtools::load_all(\"/users/mturchin/LabMisc/StephensLab/bmass\"); Data1 <- read.table(file('stdin'), header=F); for (i in 1:22) { print(i); Data2 <- Data1[Data1[,2] == i,]; Data3 <- Data2[order(-log10(Data2[,10]), decreasing=TRUE),]; Data3 <- Data3[indephits(-log10(Data2[,10]), Data2[,2], Data2[,3]) == 1,]; print(dim(Data3)); print(head(Data3)); };"
rm $HOME/data/mturchin/Data/Loh2017/body_HEIGHTz.sumstats.lt5eNeg9.bmass.GrdyClmp.txt.gz; zcat $HOME/data/mturchin/Data/Loh2017/body_HEIGHTz.sumstats.lt5eNeg9.gz | R -q -e "library(\"devtools\"); devtools::load_all(\"/users/mturchin/LabMisc/StephensLab/bmass\"); Data1 <- read.table(file('stdin'), header=F); for (i in 1:22) { print(i); Data2 <- Data1[Data1[,2] == i,]; Data3 <- Data2[order(-log10(Data2[,10]), decreasing=TRUE),]; Data3 <- Data3[indephits(-log10(Data2[,10]), Data2[,2], Data2[,3]) == 1,]; write.table(Data3, file=\"$HOME/data/mturchin/Data/Loh2017/body_HEIGHTz.sumstats.lt5eNeg9.bmass.GrdyClmp.txt\", append=TRUE, quote=FALSE, col.names=FALSE, row.names=FALSE); };"; gzip -f $HOME/data/mturchin/Data/Loh2017/body_HEIGHTz.sumstats.lt5eNeg9.bmass.GrdyClmp.txt
rm -f $HOME/data/mturchin/Data/Loh2017/body_HEIGHTz.sumstats.lt1eNeg4.bmass.GrdyClmp.txt.gz; zcat $HOME/data/mturchin/Data/Loh2017/body_HEIGHTz.sumstats.lt1eNeg4.gz | R -q -e "library(\"devtools\"); devtools::load_all(\"/users/mturchin/LabMisc/StephensLab/bmass\"); Data1 <- read.table(file('stdin'), header=F); for (i in 1:22) { print(i); Data2 <- Data1[Data1[,2] == i,]; Data3 <- Data2[order(-log10(Data2[,10]), decreasing=TRUE),]; Data3 <- Data3[indephits(-log10(Data2[,10]), Data2[,2], Data2[,3]) == 1,]; write.table(Data3, file=\"$HOME/data/mturchin/Data/Loh2017/body_HEIGHTz.sumstats.lt1eNeg4.bmass.GrdyClmp.txt\", append=TRUE, quote=FALSE, col.names=FALSE, row.names=FALSE); };"; gzip -f $HOME/data/mturchin/Data/Loh2017/body_HEIGHTz.sumstats.lt1eNeg4.bmass.GrdyClmp.txt
rm -f $HOME/data/mturchin/Data/Neale2017/Vs2/50.assoc.lt5eNeg9.bmass.GrdyClmp.tsv.gz; zcat $HOME/data/mturchin/Data/Neale2017/Vs2/50.assoc.lt5eNeg9.tsv.gz | sed 's/:/ /' | sed 's/:/ /' | R -q -e "library(\"devtools\"); devtools::load_all(\"/users/mturchin/LabMisc/StephensLab/bmass\"); Data1 <- read.table(file('stdin'), header=F); for (i in 1:22) { print(i); Data2 <- Data1[Data1[,1] == i,]; Data3 <- Data2[order(-log10(Data2[,11]), decreasing=TRUE),]; Data3 <- Data3[indephits(-log10(Data2[,11]), Data2[,1], Data2[,2]) == 1,]; write.table(Data3, file=\"$HOME/data/mturchin/Data/Neale2017/Vs2/50.assoc.lt5eNeg9.bmass.GrdyClmp.tsv\", append=TRUE, quote=FALSE, col.names=FALSE, row.names=FALSE); };"; gzip -f $HOME/data/mturchin/Data/Neale2017/Vs2/50.assoc.lt5eNeg9.bmass.GrdyClmp.tsv
rm -f $HOME/data/mturchin/Data/Neale2017/Vs2/50.assoc.lt1eNeg4.bmass.GrdyClmp.tsv.gz; zcat $HOME/data/mturchin/Data/Neale2017/Vs2/50.assoc.lt1eNeg4.tsv.gz | sed 's/:/ /' | sed 's/:/ /' | R -q -e "library(\"devtools\"); devtools::load_all(\"/users/mturchin/LabMisc/StephensLab/bmass\"); Data1 <- read.table(file('stdin'), header=F); for (i in 1:22) { print(i); Data2 <- Data1[Data1[,1] == i,]; Data3 <- Data2[order(-log10(Data2[,11]), decreasing=TRUE),]; Data3 <- Data3[indephits(-log10(Data2[,11]), Data2[,1], Data2[,2]) == 1,]; write.table(Data3, file=\"$HOME/data/mturchin/Data/Neale2017/Vs2/50.assoc.lt1eNeg4.bmass.GrdyClmp.tsv\", append=TRUE, quote=FALSE, col.names=FALSE, row.names=FALSE); };"; gzip -f $HOME/data/mturchin/Data/Neale2017/Vs2/50.assoc.lt1eNeg4.bmass.GrdyClmp.tsv
rm -f $HOME/data/mturchin/Data/GIANT2014_5/GIANT_HEIGHT_Wood_et_al_2014_publicrelease_HapMapCeuFreq.wUCSCGenomeBrowser_dbSNP130.vs1.lt5eNeg9.GrdyClmp.tx*; zcat $HOME/data/mturchin/Data/GIANT2014_5/GIANT_HEIGHT_Wood_et_al_2014_publicrelease_HapMapCeuFreq.wUCSCGenomeBrowser_dbSNP130.vs1.lt5eNeg9.txt.gz | grep -v NA | grep -v , | sed 's/_/ /g' | R -q -e "library(\"devtools\"); devtools::load_all(\"/users/mturchin/LabMisc/StephensLab/bmass\"); Data1 <- read.table(file('stdin'), header=F); for (i in 1:22) { print(i); Data2 <- Data1[Data1[,9] == i,]; Data3 <- Data2[order(-log10(Data2[,7]), decreasing=TRUE),]; Data3 <- Data3[indephits(-log10(Data2[,7]), Data2[,9], Data2[,10]) == 1,]; write.table(Data3, file=\"$HOME/data/mturchin/Data/GIANT2014_5/GIANT_HEIGHT_Wood_et_al_2014_publicrelease_HapMapCeuFreq.wUCSCGenomeBrowser_dbSNP130.vs1.lt5eNeg9.GrdyClmp.txt\", append=TRUE, quote=FALSE, col.names=FALSE, row.names=FALSE); };"; gzip -f $HOME/data/mturchin/Data/GIANT2014_5/GIANT_HEIGHT_Wood_et_al_2014_publicrelease_HapMapCeuFreq.wUCSCGenomeBrowser_dbSNP130.vs1.lt5eNeg9.GrdyClmp.txt
rm -f $HOME/data/mturchin/Data/GIANT2014_5/GIANT_HEIGHT_Wood_et_al_2014_publicrelease_HapMapCeuFreq.wUCSCGenomeBrowser_dbSNP130.vs1.lt5eNeg8.GrdyClmp.tx*; zcat $HOME/data/mturchin/Data/GIANT2014_5/GIANT_HEIGHT_Wood_et_al_2014_publicrelease_HapMapCeuFreq.wUCSCGenomeBrowser_dbSNP130.vs1.lt5eNeg8.txt.gz | grep -v NA | grep -v , | sed 's/_/ /g' | R -q -e "library(\"devtools\"); devtools::load_all(\"/users/mturchin/LabMisc/StephensLab/bmass\"); Data1 <- read.table(file('stdin'), header=F); for (i in 1:22) { print(i); Data2 <- Data1[Data1[,9] == i,]; Data3 <- Data2[order(-log10(Data2[,7]), decreasing=TRUE),]; Data3 <- Data3[indephits(-log10(Data2[,7]), Data2[,9], Data2[,10]) == 1,]; write.table(Data3, file=\"$HOME/data/mturchin/Data/GIANT2014_5/GIANT_HEIGHT_Wood_et_al_2014_publicrelease_HapMapCeuFreq.wUCSCGenomeBrowser_dbSNP130.vs1.lt5eNeg8.GrdyClmp.txt\", append=TRUE, quote=FALSE, col.names=FALSE, row.names=FALSE); };"; gzip -f $HOME/data/mturchin/Data/GIANT2014_5/GIANT_HEIGHT_Wood_et_al_2014_publicrelease_HapMapCeuFreq.wUCSCGenomeBrowser_dbSNP130.vs1.lt5eNeg8.GrdyClmp.txt
rm -f $HOME/data/mturchin/Data/GIANT2014_5/GIANT_HEIGHT_Wood_et_al_2014_publicrelease_HapMapCeuFreq.wUCSCGenomeBrowser_dbSNP130.vs1.lt1eNeg4.GrdyClmp.tx*; zcat $HOME/data/mturchin/Data/GIANT2014_5/GIANT_HEIGHT_Wood_et_al_2014_publicrelease_HapMapCeuFreq.wUCSCGenomeBrowser_dbSNP130.vs1.lt1eNeg4.txt.gz | grep -v NA | grep -v , | sed 's/_/ /g' | R -q -e "library(\"devtools\"); devtools::load_all(\"/users/mturchin/LabMisc/StephensLab/bmass\"); Data1 <- read.table(file('stdin'), header=F); for (i in 1:22) { print(i); Data2 <- Data1[Data1[,9] == i,]; Data3 <- Data2[order(-log10(Data2[,7]), decreasing=TRUE),]; Data3 <- Data3[indephits(-log10(Data2[,7]), Data2[,9], Data2[,10]) == 1,]; write.table(Data3, file=\"$HOME/data/mturchin/Data/GIANT2014_5/GIANT_HEIGHT_Wood_et_al_2014_publicrelease_HapMapCeuFreq.wUCSCGenomeBrowser_dbSNP130.vs1.lt1eNeg4.GrdyClmp.txt\", append=TRUE, quote=FALSE, col.names=FALSE, row.names=FALSE); };"; gzip -f $HOME/data/mturchin/Data/GIANT2014_5/GIANT_HEIGHT_Wood_et_al_2014_publicrelease_HapMapCeuFreq.wUCSCGenomeBrowser_dbSNP130.vs1.lt1eNeg4.GrdyClmp.txt

zcat /users/mturchin/data/mturchin/Data/Loh2017/body_HEIGHTz.sumstats.lt5eNeg9.bmass.GrdyClmp.txt.gz | awk '{ print $1 }' | gzip > /users/mturchin/data/mturchin/Data/Loh2017/body_HEIGHTz.sumstats.lt5eNeg9.bmass.GrdyClmp.rsIDs.gz 
zcat /users/mturchin/data/mturchin/Data/Loh2017/body_HEIGHTz.sumstats.lt1eNeg4.bmass.GrdyClmp.txt.gz | awk '{ print $1 }' | gzip > /users/mturchin/data/mturchin/Data/Loh2017/body_HEIGHTz.sumstats.lt1eNeg4.bmass.GrdyClmp.rsIDs.gz 
zcat /users/mturchin/data/mturchin/Data//Neale2017/Vs2/50.assoc.lt5eNeg9.bmass.GrdyClmp.tsv.gz | awk '{ print $4 }' | gzip > /users/mturchin/data/mturchin/Data//Neale2017/Vs2/50.assoc.lt5eNeg9.bmass.GrdyClmp.rsIDs.gz
zcat /users/mturchin/data/mturchin/Data//Neale2017/Vs2/50.assoc.lt1eNeg4.bmass.GrdyClmp.tsv.gz | awk '{ print $4 }' | gzip > /users/mturchin/data/mturchin/Data//Neale2017/Vs2/50.assoc.lt1eNeg4.bmass.GrdyClmp.rsIDs.gz
zcat $HOME/data/mturchin/Data/GIANT2014_5/GIANT_HEIGHT_Wood_et_al_2014_publicrelease_HapMapCeuFreq.wUCSCGenomeBrowser_dbSNP130.vs1.lt5eNeg9.GrdyClmp.txt.gz | awk '{ print $1 }' | gzip > $HOME/data/mturchin/Data/GIANT2014_5/GIANT_HEIGHT_Wood_et_al_2014_publicrelease_HapMapCeuFreq.wUCSCGenomeBrowser_dbSNP130.vs1.lt5eNeg9.GrdyClmp.rsIDs.gz 
zcat $HOME/data/mturchin/Data/GIANT2014_5/GIANT_HEIGHT_Wood_et_al_2014_publicrelease_HapMapCeuFreq.wUCSCGenomeBrowser_dbSNP130.vs1.lt5eNeg8.GrdyClmp.txt.gz | awk '{ print $1 }' | gzip > $HOME/data/mturchin/Data/GIANT2014_5/GIANT_HEIGHT_Wood_et_al_2014_publicrelease_HapMapCeuFreq.wUCSCGenomeBrowser_dbSNP130.vs1.lt5eNeg8.GrdyClmp.rsIDs.gz 
zcat $HOME/data/mturchin/Data/GIANT2014_5/GIANT_HEIGHT_Wood_et_al_2014_publicrelease_HapMapCeuFreq.wUCSCGenomeBrowser_dbSNP130.vs1.lt1eNeg4.GrdyClmp.txt.gz | awk '{ print $1 }' | gzip > $HOME/data/mturchin/Data/GIANT2014_5/GIANT_HEIGHT_Wood_et_al_2014_publicrelease_HapMapCeuFreq.wUCSCGenomeBrowser_dbSNP130.vs1.lt1eNeg4.GrdyClmp.rsIDs.gz 

#/users/mturchin/LabMisc/HirschhornLab/SohailRspnd/ukb_chrAll_v2.British.QCed.reqDrop.QCed.dropRltvs.PCAdrop.Height.Trans.ADD.assoc.linear.LohpVals.ppr.clumped.gz
#/users/mturchin/LabMisc/HirschhornLab/SohailRspnd/ukb_chrAll_v2.British.QCed.reqDrop.QCed.dropRltvs.PCAdrop.Height.Trans.ADD.assoc.linear.LohpVals.loose.clumped.gz
#/users/mturchin/LabMisc/HirschhornLab/SohailRspnd/ukb_chrAll_v2.British.QCed.reqDrop.QCed.dropRltvs.PCAdrop.Height.Trans.ADD.assoc.linear.BNpVals.ppr.clumped.gz
#/users/mturchin/LabMisc/HirschhornLab/SohailRspnd/ukb_chrAll_v2.British.QCed.reqDrop.QCed.dropRltvs.PCAdrop.Height.Trans.ADD.assoc.linear.BNpVals.loose.clumped.gz
##/users/mturchin/LabMisc/HirschhornLab/SohailRspnd/ukb_chrAll_v2.British.QCed.reqDrop.QCed.dropRltvs.PCAdrop.Height.Trans.ADD.assoc.linear.GIANTpVals.ppr.clumped.gz
##/users/mturchin/LabMisc/HirschhornLab/SohailRspnd/ukb_chrAll_v2.British.QCed.reqDrop.QCed.dropRltvs.PCAdrop.Height.Trans.ADD.assoc.linear.GIANTpVals.loose.clumped.gz
#/users/mturchin/LabMisc/HirschhornLab/20180601_Wood2014_HeightGWASSNPs.txt

mkdir /users/mturchin/data/1000G/mturchin20/Analyses
mkdir /users/mturchin/data/1000G/mturchin20/Analyses/ShiftAE
mkdir /users/mturchin/data/1000G/mturchin20/Analyses/ShiftAE/Vs1

#NOTE -- all 1000G files are on the same strand, and since pop subsets were extracted from the main source AllPops files, they should all start with the same set of SNPs at least (read: http://www.internationalgenome.org/faq/what-strand-are-variants-your-vcf-file/)

for chr1 in `echo {1..22}`; do
	echo $chr1

#	join <(join <(zcat /users/mturchin/data/1000G/subsets/CEU/AFs/CEU.chr${chr1}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.SNPs.frq.gz | awk '{ print $2 }' | sort | uniq -u) <(zcat /users/mturchin/data/1000G/subsets/CEU/AFs/CEU.chr${chr1}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.SNPs.frq.gz | awk '{ if ($6 > 125) { print $0 } } ' | awk '{ print $2 "\t" $3 "\t" $4 "\t" $5}' | sort -k 1,1)) <(join <(zcat /users/mturchin/data/1000G/subsets/GBR/AFs/GBR.chr${chr1}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.SNPs.frq.gz | awk '{ print $2 }' | sort | uniq -u) <(zcat /users/mturchin/data/1000G/subsets/GBR/AFs/GBR.chr${chr1}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.SNPs.frq.gz | awk '{ if ($6 > 125) { print $0 } } ' | awk '{ print $2 "\t" $3 "\t" $4 "\t" $5}' | sort -k 1,1)) | perl -lane 'if ($F[1] eq "0") { if ($F[4] eq "0") { if ($F[2] eq $F[5]) { print join("\t", @F); } } else { if ($F[2] eq $F[5]) { print join("\t", @F); } elsif ($F[2] eq $F[4]) { ($F[4], $F[5]) = ($F[5], $F[4]); $F[6] = 1 - $F[6]; print join("\t", @F); } else { $PH1 = 1; } } } else { if ($F[4] eq "0") { if ($F[1] eq $F[5]) { ($F[4], $F[5]) = ($F[5], $F[4]); $F[6] = 1 - $F[6]; print join("\t", @F); } elsif ($F[2] eq $F[5]) { print join("\t", @F); } else { $PH1 = 1; } } else { if (($F[1] eq $F[4]) && ($F[2] eq $F[5])) { print join("\t", @F); } elsif (($F[1] eq $F[5]) && ($F[2] eq $F[4])) { ($F[4], $F[5]) = ($F[5], $F[4]); $F[6] = 1 - $F[6]; print join("\t", @F); } else { $PH1 = 1; } } }' | \ 
#	join - <(join <(zcat /users/mturchin/data/1000G/subsets/TSI/AFs/TSI.chr${chr1}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.SNPs.frq.gz | awk '{ print $2 }' | sort | uniq -u) <(zcat /users/mturchin/data/1000G/subsets/TSI/AFs/TSI.chr${chr1}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.SNPs.frq.gz | awk '{ if ($6 > 125) { print $0 } } ' | awk '{ print $2 "\t" $3 "\t" $4 "\t" $5}' | sort -k 1,1)) | perl -lane 'if ($F[1] eq "0") { if ($F[7] eq "0") { if ($F[2] eq $F[8]) { print join("\t", @F); } } else { if ($F[2] eq $F[8]) { print join("\t", @F); } elsif ($F[2] eq $F[7]) { ($F[7], $F[8]) = ($F[8], $F[7]); $F[9] = 1 - $F[9]; print join("\t", @F); } else { $PH1 = 1; } } } else { if ($F[7] eq "0") { if ($F[1] eq $F[8]) { ($F[7], $F[8]) = ($F[8], $F[7]); $F[9] = 1 - $F[9]; print join("\t", @F); } elsif ($F[2] eq $F[8]) { print join("\t", @F); } else { $PH1 = 1; } } else { if (($F[1] eq $F[7]) && ($F[2] eq $F[8])) { print join("\t", @F); } elsif (($F[1] eq $F[8]) && ($F[2] eq $F[7])) { ($F[7], $F[8]) = ($F[8], $F[7]); $F[9] = 1 - $F[9]; print join("\t", @F); } else { $PH1 = 1; } } }' | \
#	join - <(join <(zcat /users/mturchin/data/1000G/subsets/IBS/AFs/IBS.chr${chr1}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.SNPs.frq.gz | awk '{ print $2 }' | sort | uniq -u) <(zcat /users/mturchin/data/1000G/subsets/IBS/AFs/IBS.chr${chr1}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.SNPs.frq.gz | awk '{ if ($6 > 125) { print $0 } } ' | awk '{ print $2 "\t" $3 "\t" $4 "\t" $5}' | sort -k 1,1)) | perl -lane 'if ($F[1] eq "0") { if ($F[10] eq "0") { if ($F[2] eq $F[11]) { print join("\t", @F); } } else { if ($F[2] eq $F[11]) { print join("\t", @F); } elsif ($F[2] eq $F[10]) { ($F[10], $F[11]) = ($F[11], $F[10]); $F[12] = 1 - $F[12]; print join("\t", @F); } else { $PH1 = 1; } } } else { if ($F[10] eq "0") { if ($F[1] eq $F[11]) { ($F[10], $F[11]) = ($F[11], $F[10]); $F[12] = 1 - $F[12]; print join("\t", @F); } elsif ($F[2] eq $F[11]) { print join("\t", @F); } else { $PH1 = 1; } } else { if (($F[1] eq $F[10]) && ($F[2] eq $F[11])) { print join("\t", @F); } elsif (($F[1] eq $F[11]) && ($F[2] eq $F[10])) { ($F[10], $F[11]) = ($F[11], $F[10]); $F[12] = 1 - $F[12]; print join("\t", @F); } else { $PH1 = 1; } } }' | \ 
#	join - <(join <(zcat /users/mturchin/data/1000G/subsets/YRI/AFs/YRI.chr${chr1}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.SNPs.frq.gz | awk '{ print $2 }' | sort | uniq -u) <(zcat /users/mturchin/data/1000G/subsets/YRI/AFs/YRI.chr${chr1}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.SNPs.frq.gz | awk '{ if ($6 > 125) { print $0 } } ' | awk '{ print $2 "\t" $3 "\t" $4 "\t" $5}' | sort -k 1,1)) | perl -lane 'if ($F[1] eq "0") { if ($F[13] eq "0") { if ($F[2] eq $F[14]) { print join("\t", @F); } } else { if ($F[2] eq $F[14]) { print join("\t", @F); } elsif ($F[2] eq $F[13]) { ($F[13], $F[14]) = ($F[14], $F[13]); $F[15] = 1 - $F[15]; print join("\t", @F); } else { $PH1 = 1; } } } else { if ($F[13] eq "0") { if ($F[1] eq $F[14]) { ($F[13], $F[14]) = ($F[14], $F[13]); $F[15] = 1 - $F[15]; print join("\t", @F); } elsif ($F[2] eq $F[14]) { print join("\t", @F); } else { $PH1 = 1; } } else { if (($F[1] eq $F[13]) && ($F[2] eq $F[14])) { print join("\t", @F); } elsif (($F[1] eq $F[14]) && ($F[2] eq $F[13])) { ($F[13], $F[14]) = ($F[14], $F[13]); $F[15] = 1 - $F[15]; print join("\t", @F); } else { $PH1 = 1; } } }' | \ 
#	join - <(join <(zcat /users/mturchin/data/1000G/subsets/ESN/AFs/ESN.chr${chr1}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.SNPs.frq.gz | awk '{ print $2 }' | sort | uniq -u) <(zcat /users/mturchin/data/1000G/subsets/ESN/AFs/ESN.chr${chr1}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.SNPs.frq.gz | awk '{ if ($6 > 125) { print $0 } } ' | awk '{ print $2 "\t" $3 "\t" $4 "\t" $5}' | sort -k 1,1)) | perl -lane 'if ($F[1] eq "0") { if ($F[16] eq "0") { if ($F[2] eq $F[17]) { print join("\t", @F); } } else { if ($F[2] eq $F[17]) { print join("\t", @F); } elsif ($F[2] eq $F[16]) { ($F[16], $F[17]) = ($F[17], $F[16]); $F[18] = 1 - $F[18]; print join("\t", @F); } else { $PH1 = 1; } } } else { if ($F[16] eq "0") { if ($F[1] eq $F[17]) { ($F[16], $F[17]) = ($F[17], $F[16]); $F[18] = 1 - $F[18]; print join("\t", @F); } elsif ($F[2] eq $F[17]) { print join("\t", @F); } else { $PH1 = 1; } } else { if (($F[1] eq $F[16]) && ($F[2] eq $F[17])) { print join("\t", @F); } elsif (($F[1] eq $F[17]) && ($F[2] eq $F[16])) { ($F[16], $F[17]) = ($F[17], $F[16]); $F[18] = 1 - $F[18]; print join("\t", @F); } else { $PH1 = 1; } } }' | gzip > /users/mturchin/data/1000G/mturchin20/Analyses/ShiftAE/Vs1/CEUGBRTSIESNYRIESN.chr${chr1}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.SNPs.frq.gz
#	zcat /users/mturchin/data/1000G/mturchin20/Analyses/ShiftAE/Vs1/CEUGBRTSIESNYRIESN.chr${chr1}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.SNPs.frq.gz | awk '{ if (($4 > 0) && ($7 > 0) && ($10 > 0) && ($13 > 0)) { print $0 } } ' | gzip > /users/mturchin/data/1000G/mturchin20/Analyses/ShiftAE/Vs1/CEUGBRTSIESNYRIESN.chr${chr1}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.SNPs.noEURfix.frq.gz
	zcat /users/mturchin/data/1000G/mturchin20/Analyses/ShiftAE/Vs1/CEUGBRTSIESNYRIESN.chr${chr1}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.SNPs.noEURfix.frq.gz | grep -v MAF | R -q -e "Data1 <- read.table(file('stdin'), header=F); EurMean <- apply(Data1[,c(4,7,10,13)], 1, function(x) { return(mean(x)); }); CEUDiff <- abs(Data1[,4] - EurMean); GBRDiff <- abs(Data1[,7] - EurMean); TSIDiff <- abs(Data1[,10] - EurMean); IBSDiff <- abs(Data1[,13] - EurMean); write.table(cbind(Data1[,c(1,2,3,4,7,10,13)], EurMean, CEUDiff, GBRDiff, TSIDiff, IBSDiff), quote=FALSE, col.names=FALSE, row.names=FALSE);" | grep -v \> | gzip > /users/mturchin/data/1000G/mturchin20/Analyses/ShiftAE/Vs1/CEUGBRTSIESNYRIESN.chr${chr1}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.SNPs.noEURfix.edit.wMeanInfo.frq.gz

done

mkdir /users/mturchin/data/1000G/mturchin20/Analyses/ShiftAE/Vs1/SubFiles
mkdir /users/mturchin/data/1000G/mturchin20/Analyses/ShiftAE/Vs1/SubFiles/PerAFFiles

rm -f /users/mturchin/data/1000G/mturchin20/Analyses/ShiftAE/Vs1/SubFiles/PerAFFiles/*frq*; for chr1 in `echo {1..22}`; do
	echo $chr1

	zcat /users/mturchin/data/1000G/mturchin20/Analyses/ShiftAE/Vs1/CEUGBRTSIESNYRIESN.chr${chr1}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.SNPs.noEURfix.edit.wMeanInfo.frq.gz | R -q -e "Data1 <- read.table(file('stdin'), header=F); for (i in 1:nrow(Data1)) { AncAF <- round(Data1[i,8], digits=2); write.table(Data1[i,], file=paste(\"/users/mturchin/data/1000G/mturchin20/Analyses/ShiftAE/Vs1/SubFiles/PerAFFiles/CEUGBRTSIESNYRIESN.chrAll.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.SNPs.noEURfix.edit.wMeanInfo.permPrep.AncAF_\", as.character(AncAF), \".frq\", sep=\"\"), append=TRUE, quote=FALSE, col.names=FALSE, row.names=FALSE);};"

done; gzip -f /users/mturchin/data/1000G/mturchin20/Analyses/ShiftAE/Vs1/SubFiles/PerAFFiles/*frq

mkdir /users/mturchin/data/1000G/mturchin20/Analyses/ShiftAE/Vs1/SubFiles/PerBetaFiles
mkdir /users/mturchin/data/1000G/mturchin20/Analyses/ShiftAE/Vs1/SubFiles/PerBetaFiles/Loh2017
mkdir /users/mturchin/data/1000G/mturchin20/Analyses/ShiftAE/Vs1/SubFiles/PerBetaFiles/Neale2017
mkdir /users/mturchin/data/1000G/mturchin20/Analyses/ShiftAE/Vs1/SubFiles/PerBetaFiles/GIANT2014_5

~for i in `cat <(echo "Loh2017;3169 Neale2017;9727 GIANT2014_5;27673" | perl -lane 'print join("\n", @F);') | head -n 1 | tail -n 1`; do
~	iVal1=`echo $i | perl -ane 'my @vals1 = split(/;/, $F[0]); print $vals1[0];'`; iSeed1=`echo $i | perl -ane 'my @vals1 = split(/;/, $F[0]); print $vals1[1];'`;
~	for j in `ls -lrt /users/mturchin/data/1000G/mturchin20/Analyses/ShiftAE/Vs1/SubFiles/PerAFFiles/ | awk '{ print $9 }' | head -n 10`; do
~		AFFile1="/users/mturchin/data/1000G/mturchin20/Analyses/ShiftAE/Vs1/SubFiles/PerAFFiles/${j}";
~		AFFile1NewPre=`echo "/users/mturchin/data/1000G/mturchin20/Analyses/ShiftAE/Vs1/SubFiles/PerAFFiles/${iVal1}/${j}" | sed 's/.frq.gz//g'`; AFFile1New="${AFFile1NewPre}.${iVal1}.frq.gz";
~		echo $AFFile1 $AFFile1New
~	
~		join <(zcat $AFFile1 | sort -k 1,1) <(zcat /users/mturchin/data/mturchin/Data/$iVal1/$iVal1.Height.edits.txt.gz | grep -v pval | perl -slane 'srand($iSeed2); my $IncAllele; if ($F[3] > 0) { $IncAllele = $F[1]; } elsif ($F[3] < 0) { $IncAllele =$F[2]; } elsif ($F[3] == 0) { my $rand1 = rand(); if ($rand1 > .5) { $IncAllele = $F[1]; } elsif ($rand1 <= .5) { $IncAllele = $F[2]; } else { print STDERR "Error2a -- how did this happen?"; } } else { print STDERR "Error1a -- how did this happen?"; } print $F[0], "\t", $IncAllele;' -- -iSeed2=$iSeed1 | sort -k 1,1) | \
~		perl -lane 'my $flag1 = 0; if ($F[1] eq $F[12]) { $flag1 = 1; } elsif ($F[2] eq $F[12]) { ($F[1], $F[2]) = ($F[2], $F[1]); $F[3] = 1 - $F[3]; $F[4] = 1 - $F[4]; $F[5] = 1 - $F[5]; $F[6] = 1 - $F[6]; $flag1 = 2; } else { $PH1 = 1; } print join("\t", @F), "\t", $flag1;' | \
~		R -q -e "Data1 <- read.table(file('stdin'), header=F); EurMean <- apply(Data1[,c(4:7)], 1, function(x) { return(mean(x)); }); CEUDiff <- Data1[,4] - EurMean; GBRDiff <- Data1[,5] - EurMean; TSIDiff <- Data1[,6] - EurMean; IBSDiff <- Data1[,7] - EurMean; write.table(cbind(Data1[,c(1:7)], EurMean, CEUDiff, GBRDiff, TSIDiff, IBSDiff, Data1[,ncol(Data1)-1], Data1[,ncol(Data1)]), quote=FALSE, col.names=FALSE, row.names=FALSE);" | grep -v \> | gzip > $AFFile1New
~
~	done
~done
		
for i in `cat <(echo "Loh2017;3169 Neale2017;9727 GIANT2014_5;27673" | perl -lane 'print join("\n", @F);') | head -n 1 | tail -n 1`; do
	iVal1=`echo $i | perl -ane 'my @vals1 = split(/;/, $F[0]); print $vals1[0];'`; iSeed1=`echo $i | perl -ane 'my @vals1 = split(/;/, $F[0]); print $vals1[1];'`;
	for j in `cat <(echo "Height BMI" | perl -lane 'print join("\n", @F);') | head -n 1 | tail -n 1`; do 
		for chr1 in `echo {1..22}`; do
			AFFile1="/users/mturchin/data/1000G/mturchin20/Analyses/ShiftAE/Vs1/CEUGBRTSIESNYRIESN.chr${chr1}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.SNPs.noEURfix.edit.wMeanInfo.frq.gz";
			AFFile1New="/users/mturchin/data/1000G/mturchin20/Analyses/ShiftAE/Vs1/SubFiles/PerBetaFiles/$iVal1/CEUGBRTSIESNYRIESN.chr${chr1}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.SNPs.noEURfix.edit.wMeanInfo.$iVal1.${j}.Inc.frq.gz";
			echo $AFFile1 $AFFile1New
	
			join <(zcat $AFFile1 | sort -k 1,1) <(zcat /users/mturchin/data/mturchin/Data/$iVal1/$iVal1.${j}.edits.txt.gz | grep -v pval | perl -slane 'srand($iSeed2); my $IncAllele; if ($F[3] > 0) { $IncAllele = $F[1]; } elsif ($F[3] < 0) { $IncAllele =$F[2]; } elsif ($F[3] == 0) { my $rand1 = rand(); if ($rand1 > .5) { $IncAllele = $F[1]; } elsif ($rand1 <= .5) { $IncAllele = $F[2]; } else { print STDERR "Error2a -- how did this happen?"; } } else { print STDERR "Error1a -- how did this happen?"; } print $F[0], "\t", $IncAllele;' -- -iSeed2=$iSeed1 | sort -k 1,1) | \
			perl -lane 'my $flag1 = 0; if ($F[1] eq $F[12]) { $flag1 = 1; } elsif ($F[2] eq $F[12]) { ($F[1], $F[2]) = ($F[2], $F[1]); $F[3] = 1 - $F[3]; $F[4] = 1 - $F[4]; $F[5] = 1 - $F[5]; $F[6] = 1 - $F[6]; $flag1 = 2; } else { $PH1 = 1; } print join("\t", @F), "\t", $flag1;' | \
			R -q -e "Data1 <- read.table(file('stdin'), header=F); EurMean <- apply(Data1[,c(4:7)], 1, function(x) { return(mean(x)); }); CEUDiff <- Data1[,4] - EurMean; GBRDiff <- Data1[,5] - EurMean; TSIDiff <- Data1[,6] - EurMean; IBSDiff <- Data1[,7] - EurMean; write.table(cbind(Data1[,c(1:7)], EurMean, CEUDiff, GBRDiff, TSIDiff, IBSDiff, Data1[,ncol(Data1)-1], Data1[,ncol(Data1)]), quote=FALSE, col.names=FALSE, row.names=FALSE);" | grep -v \> | gzip > $AFFile1New
		done
	done
done

mkdir /users/mturchin/data/1000G/mturchin20/Analyses/ShiftAE/Vs1/SubFiles/PerAFFiles/Loh2017
mkdir /users/mturchin/data/1000G/mturchin20/Analyses/ShiftAE/Vs1/SubFiles/PerAFFiles/Neale2017
mkdir /users/mturchin/data/1000G/mturchin20/Analyses/ShiftAE/Vs1/SubFiles/PerAFFiles/GIANT2014_5

for i in `cat <(echo "Loh2017;3169 Neale2017;9727 GIANT2014_5;27673" | perl -lane 'print join("\n", @F);') | head -n 1 | tail -n 1`; do
	iVal1=`echo $i | perl -ane 'my @vals1 = split(/;/, $F[0]); print $vals1[0];'`; iSeed1=`echo $i | perl -ane 'my @vals1 = split(/;/, $F[0]); print $vals1[1];'`;
	for j in `cat <(echo "Height BMI" | perl -lane 'print join("\n", @F);') | head -n 1 | tail -n 1`; do 
		rm -f /users/mturchin/data/1000G/mturchin20/Analyses/ShiftAE/Vs1/SubFiles/PerAFFiles/$iVal1/$j/*$j*frq*; for chr1 in `echo {1..22}`; do
		        echo $chr1
	
			if [ ! -d /users/mturchin/data/1000G/mturchin20/Analyses/ShiftAE/Vs1/SubFiles/PerAFFiles/$iVal1/$j ] ; then
				mkdir /users/mturchin/data/1000G/mturchin20/Analyses/ShiftAE/Vs1/SubFiles/PerAFFiles/$iVal1/$j
			fi
	
		        zcat /users/mturchin/data/1000G/mturchin20/Analyses/ShiftAE/Vs1/SubFiles/PerBetaFiles/$iVal1/CEUGBRTSIESNYRIESN.chr${chr1}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.SNPs.noEURfix.edit.wMeanInfo.$iVal1.${j}.Inc.frq.gz | R -q -e "Data1 <- read.table(file('stdin'), header=F); for (i in 1:nrow(Data1)) { AncAF <- round(Data1[i,8], digits=2); write.table(Data1[i,], file=paste(\"/users/mturchin/data/1000G/mturchin20/Analyses/ShiftAE/Vs1/SubFiles/PerAFFiles/$iVal1/$j/CEUGBRTSIESNYRIESN.chrAll.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.SNPs.noEURfix.edit.wMeanInfo.$iVal1.${j}.Inc.permPrep.AncAF_\", as.character(AncAF), \".frq\", sep=\"\"), append=TRUE, quote=FALSE, col.names=FALSE, row.names=FALSE);};"
		
		done; gzip -f /users/mturchin/data/1000G/mturchin20/Analyses/ShiftAE/Vs1/SubFiles/PerAFFiles/$iVal1/$j/*$j*frq
	done;
done;

##mkdir /users/mturchin/LabMisc/RamachandranLab/InterPath/ShiftAE/Data
##mkdir /users/mturchin/LabMisc/RamachandranLab/InterPath/ShiftAE/Data/GWASsnps
mkdir /users/mturchin/LabMisc/RamachandranLab/ShiftAE/Data
mkdir /users/mturchin/LabMisc/RamachandranLab/ShiftAE/Data/GWASsnps

ln -s /users/mturchin/data/mturchin/Data/Loh2017/body_HEIGHTz.sumstats.lt5eNeg9.bmass.GrdyClmp.rsIDs.gz /users/mturchin/LabMisc/RamachandranLab/ShiftAE/Data/GWASsnps/Loh2017.Height.lt5eNeg9.bmass.GrdyClmp.rsIDs.gz
ln -s /users/mturchin/data/mturchin/Data/Loh2017/body_HEIGHTz.sumstats.lt1eNeg4.bmass.GrdyClmp.rsIDs.gz /users/mturchin/LabMisc/RamachandranLab/ShiftAE/Data/GWASsnps/Loh2017.Height.lt1eNeg4.bmass.GrdyClmp.rsIDs.gz
ln -s /users/mturchin/data/mturchin/Data/Neale2017/Vs2/50.assoc.lt5eNeg9.bmass.GrdyClmp.rsIDs.gz /users/mturchin/LabMisc/RamachandranLab/ShiftAE/Data/GWASsnps/Neale2017.Height.lt5eNeg9.bmass.GrdyClmp.rsIDs.gz
ln -s /users/mturchin/data/mturchin/Data/Neale2017/Vs2/50.assoc.lt1eNeg4.bmass.GrdyClmp.rsIDs.gz /users/mturchin/LabMisc/RamachandranLab/ShiftAE/Data/GWASsnps/Neale2017.Height.lt1eNeg4.bmass.GrdyClmp.rsIDs.gz
ln -s $HOME/data/mturchin/Data/GIANT2014_5/GIANT_HEIGHT_Wood_et_al_2014_publicrelease_HapMapCeuFreq.wUCSCGenomeBrowser_dbSNP130.vs1.lt5eNeg9.GrdyClmp.rsIDs.gz /users/mturchin/LabMisc/RamachandranLab/ShiftAE/Data/GWASsnps/GIANT2014_5.Height.lt5eNeg9.bmass.GrdyClmp.rsIDs.gz
ln -s $HOME/data/mturchin/Data/GIANT2014_5/GIANT_HEIGHT_Wood_et_al_2014_publicrelease_HapMapCeuFreq.wUCSCGenomeBrowser_dbSNP130.vs1.lt5eNeg8.GrdyClmp.rsIDs.gz /users/mturchin/LabMisc/RamachandranLab/ShiftAE/Data/GWASsnps/GIANT2014_5.Height.lt5eNeg8.bmass.GrdyClmp.rsIDs.gz
ln -s $HOME/data/mturchin/Data/GIANT2014_5/GIANT_HEIGHT_Wood_et_al_2014_publicrelease_HapMapCeuFreq.wUCSCGenomeBrowser_dbSNP130.vs1.lt1eNeg4.GrdyClmp.rsIDs.gz /users/mturchin/LabMisc/RamachandranLab/ShiftAE/Data/GWASsnps/GIANT2014_5.Height.lt1eNeg4.bmass.GrdyClmp.rsIDs.gz

for i in `cat <(echo "Loh2017 Neale2017 GIANT2014_5" | perl -lane 'print join("\n", @F);') | head -n 3 | tail -n 3`; do
	for j in `cat <(echo "Height BMI" | perl -lane 'print join("\n", @F);') | head -n 1 | tail -n 1`; do
		for k in `cat <(echo "lt5eNeg9 lt1eNeg4 lt5eNeg8" | perl -lane 'print join("\n", @F);') | head -n 2 | tail -n 2`; do
			echo $i $j $k

##			join <(zcat /users/mturchin/LabMisc/RamachandranLab/InterPath/ShiftAE/Data/GWASsnps/$i.$j.$k.bmass.GrdyClmp.rsIDs.gz | sort) <(zcat /users/mturchin/data/1000G/mturchin20/Analyses/ShiftAE/Vs1/CEUGBRTSIESNYRIESN.chr*.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.SNPs.noEURfix.edit.wMeanInfo.frq.gz | sort -k 1,1) | gzip > /users/mturchin/LabMisc/RamachandranLab/InterPath/ShiftAE/Data/GWASsnps/$i.$j.$k.bmass.GrdyClmp.rsIDs.w1000GInfo.txt.gz
##			zcat /users/mturchin/LabMisc/RamachandranLab/InterPath/ShiftAE/Data/GWASsnps/$i.$j.$k.bmass.GrdyClmp.rsIDs.w1000GInfo.txt.gz | R -q -e "Data1 <- read.table(file('stdin'), header=F); write.table(matrix(c(sum(abs(mean(Data1[,ncol(Data1)-3]) - Data1[,ncol(Data1)-3])), sum(abs(mean(Data1[,ncol(Data1)-2]) - Data1[,ncol(Data1)-2])), sum(abs(mean(Data1[,ncol(Data1)-1]) - Data1[,ncol(Data1)-1])), sum(abs(mean(Data1[,ncol(Data1)]) - Data1[,ncol(Data1)]))), nrow=1), quote=FALSE, col.names=FALSE, row.names=FALSE);"	
			join <(zcat /users/mturchin/LabMisc/RamachandranLab/ShiftAE/Data/GWASsnps/$i.$j.$k.bmass.GrdyClmp.rsIDs.gz | sort) <(zcat /users/mturchin/data/1000G/mturchin20/Analyses/ShiftAE/Vs1/SubFiles/PerBetaFiles/$i/CEUGBRTSIESNYRIESN.chr*.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.SNPs.noEURfix.edit.wMeanInfo.${i}.${j}.Inc.frq.gz | sort -k 1,1) | gzip > /users/mturchin/LabMisc/RamachandranLab/ShiftAE/Data/GWASsnps/$i.$j.$k.bmass.GrdyClmp.rsIDs.w1000GInfo.${j}.Inc.txt.gz

 /users/mturchin/data/1000G/mturchin20/Analyses/ShiftAE/Vs1/CEUGBRTSIESNYRIESN.chr*.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.SNPs.noEURfix.edit.wMeanInfo.frq.gz | sort -k 1,1) | gzip > /users/mturchin/LabMisc/RamachandranLab/InterPath/ShiftAE/Data/GWASsnps/$i.$j.$k.bmass.GrdyClmp.rsIDs.w1000GInfo.txt.gz

/users/mturchin/data/1000G/mturchin20/Analyses/ShiftAE/Vs1/SubFiles/PerBetaFiles/$iVal1/CEUGBRTSIESNYRIESN.chr${chr1}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.SNPs.noEURfix.edit.wMeanInfo.$iVal1.${j}.Inc.frq.gz

		done
	done
done

mkdir /users/mturchin/data/1000G/mturchin20/Analyses/ShiftAE/Vs1/SubFiles/PermFiles

#From: https://stackoverflow.com/questions/14904983/how-do-i-check-the-existence-of-a-local-file
for i in `cat <(echo "Loh2017;3169 Neale2017;9727 GIANT2014_5;27673" | perl -lane 'print join("\n", @F);') | head -n 1 | tail -n 1`; do
	iVal1=`echo $i | perl -ane 'my @vals1 = split(/;/, $F[0]); print $vals1[0];'`; iSeed1=`echo $i | perl -ane 'my @vals1 = split(/;/, $F[0]); print $vals1[1];'`;
	for j in `cat <(echo "Height;5638 BMI;1456" | perl -lane 'print join("\n", @F);') | head -n 1 | tail -n 1`; do
		jVal1=`echo $j | perl -ane 'my @vals1 = split(/;/, $F[0]); print $vals1[0];'`; jSeed1=`echo $j | perl -ane 'my @vals1 = split(/;/, $F[0]); print $vals1[1];'`;
		for k in `cat <(echo "lt5eNeg9;2759 lt1eNeg4;78364 lt5eNeg8;3869" | perl -lane 'print join("\n", @F);') | head -n 1 | tail -n 1`; do
			kVal1=`echo $k | perl -ane 'my @vals1 = split(/;/, $F[0]); print $vals1[0];'`; kSeed1=`echo $k | perl -ane 'my @vals1 = split(/;/, $F[0]); print $vals1[1];'`;
			TotalSeed1=$((iSeed1+jSeed1+kSeed1));
			echo $i $j $k $TotalSeed1

			if [ ! -d /users/mturchin/data/1000G/mturchin20/Analyses/ShiftAE/Vs1/SubFiles/PermFiles/$iVal1/$jVal1/$kVal1 ]; then
				mkdir /users/mturchin/data/1000G/mturchin20/Analyses/ShiftAE/Vs1/SubFiles/PermFiles/$iVal1
				mkdir /users/mturchin/data/1000G/mturchin20/Analyses/ShiftAE/Vs1/SubFiles/PermFiles/$iVal1/$jVal1
				mkdir /users/mturchin/data/1000G/mturchin20/Analyses/ShiftAE/Vs1/SubFiles/PermFiles/$iVal1/$jVal1/$kVal1
			fi

			rm -f /users/mturchin/data/1000G/mturchin20/Analyses/ShiftAE/Vs1/SubFiles/PermFiles/$iVal1/$jVal1/$kVal1/*
			zcat /users/mturchin/LabMisc/RamachandranLab/ShiftAE/Data/GWASsnps/$iVal1.$jVal1.$kVal1.bmass.GrdyClmp.rsIDs.w1000GInfo.txt.gz | R -q -e "Data1 <- read.table(file('stdin'), header=F); \
			for (i in 1:10) { \
				AncAF <- round(as.numeric(as.character(Data1[i,8])), digits=2); Data2 <- c(); File3Count <- c(); for (j in c(-.01,0,.01)) { \
		    Filename3 <- paste(\"/users/mturchin/data/1000G/mturchin20/Analyses/ShiftAE/Vs1/SubFiles/PerAFFiles/$iVal1/$jVal1/CEUGBRTSIESNYRIESN.chrAll.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.SNPs.noEURfix.edit.wMeanInfo.$iVal1.${jVal1}.Inc.permPrep.AncAF_\", as.character(AncAF + j), \".frq.gz\", sep=\"\"); if (file.exists(Filename3)) { Data3 <- read.table(Filename3, header=F); Data2 <- rbind(Data2, Data3); File3Count <- c(File3Count, j); }; }; \
				cat(c(i, as.numeric(as.character($TotalSeed1))+i, AncAF, paste(File3Count, collapse=\",\"), dim(Data2)), \"\n\", sep=\"\t\"); \ 
				set.seed(as.numeric(as.character($TotalSeed1))+i); RowVals1 <- sample(1:nrow(Data2)); for (k in 1:20) { \
						write.table(Data2[RowVals1[k],], file=paste(\"/users/mturchin/data/1000G/mturchin20/Analyses/ShiftAE/Vs1/SubFiles/PermFiles/$iVal1/$jVal1/$kVal1/CEUGBRTSIESNYRIESN.chrAll.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.SNPs.noEURfix.edit.wMeanInfo.$iVal1.${jVal1}.Inc.$kVal1.perm\", as.character(k), \".txt\", sep=\"\"), append=TRUE, quote=FALSE, row.names=FALSE, col.names=FALSE); \
				}; \
			};" | grep -v ^\> > /users/mturchin/data/1000G/mturchin20/Analyses/ShiftAE/Vs1/SubFiles/PermFiles/$iVal1/$jVal1/$kVal1/CEUGBRTSIESNYRIESN.chrAll.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.SNPs.noEURfix.edit.wMeanInfo.$iVal1.${jVal1}.Inc.$kVal1.permAll.Summary.txt
			gzip -f /users/mturchin/data/1000G/mturchin20/Analyses/ShiftAE/Vs1/SubFiles/PermFiles/$iVal1/$jVal1/$kVal1/*

		done
	done
done

##mkdir /users/mturchin/data/1000G/mturchin20/Analyses/ShiftAE/Vs1/Analyses
##mkdir /users/mturchin/data/1000G/mturchin20/Analyses/ShiftAE/Vs1/Analyses/Perms
mkdir /users/mturchin/LabMisc/RamachandranLab/InterPath/ShiftAE/Vs1/Analyses/Perms

for i in `cat <(echo "Loh2017;3169 Neale2017;9727 GIANT2014_5;27673" | perl -lane 'print join("\n", @F);') | head -n 3 | tail -n 3`; do
        iVal1=`echo $i | perl -ane 'my @vals1 = split(/;/, $F[0]); print $vals1[0];'`; iSeed1=`echo $i | perl -ane 'my @vals1 = split(/;/, $F[0]); print $vals1[1];'`;
        for j in `cat <(echo "Height;5638" | perl -lane 'print join("\n", @F);') | head -n 1 | tail -n 1`; do
                jVal1=`echo $j | perl -ane 'my @vals1 = split(/;/, $F[0]); print $vals1[0];'`; jSeed1=`echo $j | perl -ane 'my @vals1 = split(/;/, $F[0]); print $vals1[1];'`;
                for k in `cat <(echo "lt5eNeg9;2759 lt1eNeg4;78364 lt5eNeg8;3869" | perl -lane 'print join("\n", @F);') | head -n 2 | tail -n 2`; do
                        kVal1=`echo $k | perl -ane 'my @vals1 = split(/;/, $F[0]); print $vals1[0];'`; kSeed1=`echo $k | perl -ane 'my @vals1 = split(/;/, $F[0]); print $vals1[1];'`;
                        TotalSeed1=$((iSeed1+jSeed1+kSeed1));
                        echo $i $j $k $TotalSeed1

			rm -f /users/mturchin/LabMisc/RamachandranLab/InterPath/ShiftAE/Vs1/Analyses/Perms/$iVal1.$jVal1.$kVal1.CEUGBRTSIESNYRIESN.chrAll.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.SNPs.noEURfix.edit.wMeanInfo.permAll.Results.txt; for l in {1..20}; do
				zcat /users/mturchin/data/1000G/mturchin20/Analyses/ShiftAE/Vs1/SubFiles/PermFiles/$iVal1.$jVal1.$kVal1/CEUGBRTSIESNYRIESN.chrAll.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.SNPs.noEURfix.edit.wMeanInfo.perm${l}.txt.gz | R -q -e "Data1 <- read.table(file('stdin'), header=F); write.table(matrix(c(sum(abs(mean(Data1[,ncol(Data1)-3]) - Data1[,ncol(Data1)-3])), sum(abs(mean(Data1[,ncol(Data1)-2]) - Data1[,ncol(Data1)-2])), sum(abs(mean(Data1[,ncol(Data1)-1]) - Data1[,ncol(Data1)-1])), sum(abs(mean(Data1[,ncol(Data1)]) - Data1[,ncol(Data1)]))), nrow=1), quote=FALSE, col.names=FALSE, row.names=FALSE);" | grep -v ^\> >> /users/mturchin/LabMisc/RamachandranLab/InterPath/ShiftAE/Vs1/Analyses/Perms/$iVal1.$jVal1.$kVal1.CEUGBRTSIESNYRIESN.chrAll.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.SNPs.noEURfix.edit.wMeanInfo.permAll.Results.txt
			done 

                done
        done
done

#				...Data1[,ncol(Data1)-3]), sd(Data1[,ncol(Data1)-2]), sd(Data1[,ncol(Data1)-1]), sd(Data1[,ncol(Data1)])), nrow=1), quote=FALSE, col.names=FALSE, row.names=FALSE);" | grep -v ^\> >> /users/mturchin/LabMisc/RamachandranLab/InterPath/ShiftAE/Vs1/Analyses/Perms/$iVal1.$jVal1.$kVal1.CEUGBRTSIESNYRIESN.chrAll.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.SNPs.noEURfix.edit.wMeanInfo.permAll.Results.txt 

for i in `cat <(echo "Loh2017;3169 Neale2017;9727 GIANT2014_5;27673" | perl -lane 'print join("\n", @F);') | head -n 3 | tail -n 3`; do
        iVal1=`echo $i | perl -ane 'my @vals1 = split(/;/, $F[0]); print $vals1[0];'`; iSeed1=`echo $i | perl -ane 'my @vals1 = split(/;/, $F[0]); print $vals1[1];'`;
        for j in `cat <(echo "Height;5638" | perl -lane 'print join("\n", @F);') | head -n 1 | tail -n 1`; do
                jVal1=`echo $j | perl -ane 'my @vals1 = split(/;/, $F[0]); print $vals1[0];'`; jSeed1=`echo $j | perl -ane 'my @vals1 = split(/;/, $F[0]); print $vals1[1];'`;
                for k in `cat <(echo "lt5eNeg9;2759 lt1eNeg4;78364 lt5eNeg8;3869" | perl -lane 'print join("\n", @F);') | head -n 2 | tail -n 2`; do
                        kVal1=`echo $k | perl -ane 'my @vals1 = split(/;/, $F[0]); print $vals1[0];'`; kSeed1=`echo $k | perl -ane 'my @vals1 = split(/;/, $F[0]); print $vals1[1];'`;
                        TotalSeed1=$((iSeed1+jSeed1+kSeed1));
                        echo $i $j $k $TotalSeed1

			cat /users/mturchin/LabMisc/RamachandranLab/InterPath/ShiftAE/Vs1/Analyses/Perms/$iVal1.$jVal1.$kVal1.CEUGBRTSIESNYRIESN.chrAll.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.SNPs.noEURfix.edit.wMeanInfo.permAll.Results.txt | head -n 10

                done
        done
done

















~~~
#20190324
[  mturchin@node1303  ~/LabMisc/RamachandranLab/InterPath/CAFE]$join <(zcat /users/mturchin/data/1000G/subsets/CEU/AFs/CEU.chr${chr1}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.SNPs.frq.gz | awk '{ print $2 "\t" $3 "\t" $4 "\t" $5}' | sort -k 1,1) <(zcat /users/mturchin/data/1000G/subsets/GBR/AFs/GBR.chr${chr1}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.SNPs.frq.gz | awk '{ print $2 "\t" $3 "\t" $4 "\t" $5}' | sort -k 1,1) | perl -lane 'if ($F[1] eq "0") { if ($F[4] eq "0") { if ($F[2] eq $F[5]) { print join("\t", @F); } } else { if ($F[2] eq $F[5]) { print join("\t", @F); } elsif ($F[2] eq $F[4]) { ($F[4], $F[5]) = ($F[5], $F[4]); $F[6] = 1 - $F[6]; print join("\t", @F); } else { $PH1 = 1; } } } else { if ($F[4] eq "0") { if ($F[1] eq $F[5]) { ($F[4], $F[5]) = ($F[5], $F[4]); $F[6] = 1 - $F[6]; print join("\t", @F); } elsif ($F[2] eq $F[5]) { print join("\t", @F); } else { $PH1 = 1; } } else { if (($F[1] eq $F[4]) && ($F[2] eq $F[5])) { print join("\t", @F); } elsif (($F[1] eq $F[5]) && ($F[2] eq $F[4])) { ($F[4], $F[5]) = ($F[5], $F[4]); $F[6] = 1 - $F[6]; print join("\t", @F); } else { $PH1 = 1; } } }' | head -n 10
SNP     A1      A2      MAF     A1      A2      MAF
rs1000033       G       T       0.2071  G       T       0.2033
rs1000050       C       T       0.1364  C       T       0.1429
rs1000070       T       C       0.2828  T       C       0.3626
rs1000072       0       G       0       0       G       0
rs1000073       A       G       0.3283  A       G       0.4396
rs1000075       T       C       0.3586  T       C       0.2473
rs1000085       C       G       0.2121  C       G       0.2308
rs1000127       C       T       0.3081  C       T       0.2967
rs1000184       C       G       0.2879  C       G       0.2253
[  mturchin@node1303  ~/LabMisc/RamachandranLab/InterPath/CAFE]$join <(zcat /users/mturchin/data/1000G/subsets/CEU/AFs/CEU.chr${chr1}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.SNPs.frq.gz | awk '{ print $2 "\t" $3 "\t" $4 "\t" $5}' | sort -k 1,1) <(zcat /users/mturchin/data/1000G/subsets/GBR/AFs/GBR.chr${chr1}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.SNPs.frq.gz | awk '{ print $2 "\t" $3 "\t" $4 "\t" $5}' | sort -k 1,1) | perl -lane 'if ($F[1] eq "0") { if ($F[4] eq "0") { if ($F[2] eq $F[5]) { print join("\t", @F); } } else { if ($F[2] eq $F[5]) { print join("\t", @F); } elsif ($F[2] eq $F[4]) { ($F[4], $F[5]) = ($F[5], $F[4]); $F[6] = 1 - $F[6]; print join("\t", @F); } else { $PH1 = 1; } } } else { if ($F[4] eq "0") { if ($F[1] eq $F[5]) { ($F[4], $F[5]) = ($F[5], $F[4]); $F[6] = 1 - $F[6]; print join("\t", @F); } elsif ($F[2] eq $F[5]) { print join("\t", @F); } else { $PH1 = 1; } } else { if (($F[1] eq $F[4]) && ($F[2] eq $F[5])) { print join("\t", @F); } elsif (($F[1] eq $F[5]) && ($F[2] eq $F[4])) { ($F[4], $F[5]) = ($F[5], $F[4]); $F[6] = 1 - $F[6]; print join("\t", @F); } else { $PH1 = 1; } } }' | wc        
6196192 43373344 158292573
[  mturchin@node1303  ~/LabMisc/RamachandranLab/InterPath/CAFE]$zcat /users/mturchin/data/1000G/subsets/CEU/AFs/CEU.chr${chr1}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.SNPs.frq.gz | wc
6196152 37176912 316003752
[  mturchin@node1303  ~/LabMisc/RamachandranLab/InterPath/CAFE]$zcat /users/mturchin/data/1000G/subsets/GBR/AFs/GBR.chr${chr1}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.SNPs.frq.gz | wc
6196152 37176912 316003752
#20190616
[  mturchin@node1137  ~/LabMisc/RamachandranLab/InterPath/CAFE]$       zcat /users/mturchin/data/1000G/mturchin20/Analyses/CAFE/Vs1/CEUGBRTSIESNYRIESN.chr${chr1}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.SNPs.noEURfix.frq.gz | R -q -e "Data1 <- read.table(file('stdin'), header=T); EurMean <- apply(Data1[,c(4,7,10,13)], 1, function(x) { return(mean(x)); }); CEUDiff <- abs(Data1[,4] - EurMean); GBRDiff <- abs(Data1[,7] - EurMean); TSIDiff <- abs(Data1[,10] - EurMean); IBSDiff <- abs(Data1[,13] - EurMean); write.table(cbind(Data1[,c(1,2,3,4,7,10,13)], EurMean, CEUDiff, GBRDiff, TSIDiff, IBSDiff), quote=FALSE, col.names=FALSE, row.names=FALSE);" | grep -v \> | gzip > /users/mturchin/data/1000G/mturchin20/Analyses/CAFE/Vs1/CEUGBRTSIESNYRIESN.chr${chr1}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.SNPs.noEURfix.edit.wMeanInfo.frq.gz
[  mturchin@node1137  ~/LabMisc/RamachandranLab/InterPath/CAFE]$zcat /users/mturchin/data/1000G/mturchin20/Analyses/CAFE/Vs1/CEUGBRTSIESNYRIESN.chr${chr1}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.SNPs.noEURfix.edit.wMeanInfo.frq.gz | head -n 10
rs1000033 G T 0.2071 0.2033 0.1682 0.1589 0.184375 0.022725 0.018925 0.016175 0.025475
rs1000050 C T 0.1364 0.1429 0.1449 0.215 0.1598 0.0234 0.0169 0.0149 0.0552
rs1000070 T C 0.2828 0.3626 0.3224 0.3178 0.3214 0.0386 0.0412 0.001 0.00359999999999999
rs1000073 A G 0.3283 0.4396 0.4953 0.4019 0.416275 0.087975 0.023325 0.079025 0.014375
rs1000075 T C 0.3586 0.2473 0.3458 0.3131 0.3162 0.0424 0.0689 0.0296 0.00309999999999999
rs1000085 C G 0.2121 0.2308 0.1869 0.2196 0.21235 0.00025 0.01845 0.02545 0.00724999999999998
rs1000127 C T 0.3081 0.2967 0.2664 0.2991 0.292575 0.015525 0.00412499999999999 0.026175 0.00652499999999995
rs1000184 C G 0.2879 0.2253 0.1776 0.229 0.22995 0.05795 0.00464999999999999 0.05235 0.000949999999999979
rs1000282 C T 0.197 0.2363 0.1916 0.2009 0.20645 0.00944999999999999 0.02985 0.01485 0.00555
rs1000283 A G 0.197 0.2363 0.1916 0.2009 0.20645 0.00944999999999999 0.02985 0.01485 0.00555
[  mturchin@node1137  ~/LabMisc/RamachandranLab/InterPath/CAFE]$zcat /users/mturchin/data/1000G/mturchin20/Analyses/CAFE/Vs1/CEUGBRTSIESNYRIESN.chr${chr1}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.SNPs.noEURfix.edit.wMeanInfo.frq.gz | wc        
 668049 8016588 64007801
[  mturchin@node1137  ~/LabMisc/RamachandranLab/InterPath/CAFE]$zcat /users/mturchin/data/1000G/mturchin20/Analyses/CAFE/Vs1/CEUGBRTSIESNYRIESN.chr${chr1}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.SNPs.noEURfix.frq.gz | wc
 668050 12692950 50935772
[  mturchin@node1137  ~/LabMisc/RamachandranLab/InterPath/CAFE]$zcat /users/mturchin/data/1000G/mturchin20/Analyses/CAFE/Vs1/CEUGBRTSIESNYRIESN.chr${chr1}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.SNPs.noEURfix.frq.gz | head -n 10
SNP     A1      A2      MAF     A1      A2      MAF     A1      A2      MAF     A1      A2      MAF     A1      A2      MAF     A1      A2      MAF
rs1000033       G       T       0.2071  G       T       0.2033  G       T       0.1682  G       T       0.1589  G       T       0.3843  G       T       0.4293
rs1000050       C       T       0.1364  C       T       0.1429  C       T       0.1449  C       T       0.215   C       T       0.625   C       T       0.6616
rs1000070       T       C       0.2828  T       C       0.3626  T       C       0.3224  T       C       0.3178  T       C       0.6528  T       C       0.6313
rs1000073       A       G       0.3283  A       G       0.4396  A       G       0.4953  A       G       0.4019  A       G       0.5741  A       G       0.6061
rs1000075       T       C       0.3586  T       C       0.2473  T       C       0.3458  T       C       0.3131  T       C       0.2685  T       C       0.3232
rs1000085       C       G       0.2121  C       G       0.2308  C       G       0.1869  C       G       0.2196  C       G       0.01389 C       G       0.07576
rs1000127       C       T       0.3081  C       T       0.2967  C       T       0.2664  C       T       0.2991  C       T       0.1852  C       T       0.1566
rs1000184       C       G       0.2879  C       G       0.2253  C       G       0.1776  C       G       0.229   C       G       0.00463 0       G       0
rs1000282       C       T       0.197   C       T       0.2363  C       T       0.1916  C       T       0.2009  C       T       0.1667  C       T       0.1364
[  mturchin@node1162  ~/data/mturchin/Loh2017]$zcat /users/mturchin/data/1000G/mturchin20/Analyses/CAFE/Vs1/SubFiles/PermFiles/Loh2017.Height.lt5eNeg9/CEUGBRTSIESNYRIESN.chrAll.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.SNPs.noEURfix.edit.wMeanInfo.perm1.txt.gz
rs12673771 C T 0.4949 0.5659 0.528 0.5327 0.530375 0.035475 0.0355249999999999 0.00237500000000002 0.00232499999999991
rs75706943 T C 0.1566 0.1264 0.1729 0.1355 0.14785 0.00874999999999998 0.02145 0.02505 0.01235
rs111252186 G A 0.06566 0.05495 0.04206 0.04206 0.0511825 0.0144775 0.0037675 0.0091225 0.0091225
rs62511744 G T 0.1414 0.1978 0.1589 0.1682 0.166575 0.025175 0.031225 0.00767499999999999 0.00162499999999999
rs3904796 T C 0.2727 0.3571 0.271 0.3411 0.310475 0.037775 0.046625 0.039475 0.030625
rs7691417 T C 0.3939 0.511 0.4626 0.4299 0.44935 0.0554500000000001 0.06165 0.01325 0.01945
rs2757579 G T 0.3788 0.2967 0.3598 0.3318 0.341775 0.037025 0.045075 0.018025 0.00997500000000001
rs1731433 G A 0.1717 0.1154 0.1729 0.1869 0.161725 0.00997499999999998 0.046325 0.011175 0.025175
rs139664127 A T 0.02525 0.01648 1 1 0.5104325 0.4851825 0.4939525 0.4895675 0.4895675
rs77942014 A G 0.1061 0.0989 0.1028 0.06075 0.0921375 0.0139625 0.0067625 0.0106625 0.0313875
#Change seed value to a constant
[  mturchin@node1162  ~/data/mturchin/Loh2017]$zcat /users/mturchin/data/1000G/mturchin20/Analyses/CAFE/Vs1/SubFiles/PermFiles/Loh2017.Height.lt5eNeg9/CEUGBRTSIESNYRIESN.chrAll.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.SNPs.noEURfix.edit.wMeanInfo.perm1.txt.gz
rs1117507 G A 0.4848 0.4945 0.6028 0.5187 0.5252 0.0404 0.0307 0.0776 0.00649999999999995
rs72963507 A G 0.1212 0.1319 0.1729 0.1542 0.14505 0.02385 0.01315 0.02785 0.00914999999999999
rs75479897 A G 0.05051 0.04945 0.04206 0.04206 0.04602 0.00449 0.00343 0.00396 0.00396
rs6822200 C T 0.1212 0.1978 0.1776 0.1916 0.17205 0.05085 0.02575 0.00555 0.01955
rs6446530 C T 0.303 0.3462 0.3178 0.2664 0.30835 0.00535000000000002 0.03785 0.00945000000000001 0.04195
rs7654800 G C 0.4192 0.4176 0.472 0.472 0.4452 0.026 0.0276 0.0268 0.0268
rs6847767 C G 0.3333 0.3462 0.3364 0.3458 0.340425 0.00712499999999999 0.00577500000000003 0.004025 0.00537500000000002
rs72648749 T C 0.1515 0.1648 0.1542 0.1495 0.155 0.0035 0.0098 0.000799999999999995 0.0055
rs6830604 C T 0.4697 0.4835 0.5654 0.5187 0.509325 0.039625 0.025825 0.056075 0.00937500000000002
rs76269601 T C 0.07071 0.09341 0.08411 0.1075 0.0889325 0.0182225 0.00447750000000001 0.00482249999999999 0.0185675
#Change seed value back to original code
[  mturchin@node1162  ~/data/mturchin/Loh2017]$zcat /users/mturchin/data/1000G/mturchin20/Analyses/CAFE/Vs1/SubFiles/PermFiles/Loh2017.Height.lt5eNeg9/CEUGBRTSIESNYRIESN.chrAll.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.SNPs.noEURfix.edit.wMeanInfo.perm1.txt.gz 
rs12673771 C T 0.4949 0.5659 0.528 0.5327 0.530375 0.035475 0.0355249999999999 0.00237500000000002 0.00232499999999991
rs75706943 T C 0.1566 0.1264 0.1729 0.1355 0.14785 0.00874999999999998 0.02145 0.02505 0.01235
rs111252186 G A 0.06566 0.05495 0.04206 0.04206 0.0511825 0.0144775 0.0037675 0.0091225 0.0091225
rs62511744 G T 0.1414 0.1978 0.1589 0.1682 0.166575 0.025175 0.031225 0.00767499999999999 0.00162499999999999
rs3904796 T C 0.2727 0.3571 0.271 0.3411 0.310475 0.037775 0.046625 0.039475 0.030625
rs7691417 T C 0.3939 0.511 0.4626 0.4299 0.44935 0.0554500000000001 0.06165 0.01325 0.01945
rs2757579 G T 0.3788 0.2967 0.3598 0.3318 0.341775 0.037025 0.045075 0.018025 0.00997500000000001
rs1731433 G A 0.1717 0.1154 0.1729 0.1869 0.161725 0.00997499999999998 0.046325 0.011175 0.025175
rs139664127 A T 0.02525 0.01648 1 1 0.5104325 0.4851825 0.4939525 0.4895675 0.4895675
rs77942014 A G 0.1061 0.0989 0.1028 0.06075 0.0921375 0.0139625 0.0067625 0.0106625 0.0313875
#20200111
[  mturchin@node1149  ~]$zcat $HOME/data/mturchin/Neale2017/Vs2/50.assoc.tsv.gz | sed 's/:/ /g' | head -n 10
variant rsid    nCompleteSamples        AC      ytx     beta    se      tstat   pval
5 43888254 C T  rs13184706      336474  1.23213e+04     3.83186e+02     -1.33660e-02    6.48300e-03     -2.06170e+00    3.92375e-02
5 43888493 C T  rs58824264      336474  2.42483e+03     3.70381e+01     -1.88438e-02    1.44299e-02     -1.30589e+00    1.91592e-01
5 43888556 T C  rs72762387      336474  1.64428e+04     6.65756e+02     -9.45691e-03    5.62698e-03     -1.68064e+00    9.28345e-02
5 43888648 C T  rs115032754     336474  1.35047e+04     5.66843e+02     1.06178e-03     6.29484e-03     1.68674e-01     8.66053e-01
5 43888690 C G  rs147555725     336474  1.24755e+03     4.51586e+01     -4.19957e-03    2.06522e-02     -2.03347e-01    8.38864e-01
5 43888838 G C  rs13185925      336474  2.33424e+04     8.95665e+02     -8.69909e-03    4.74596e-03     -1.83294e+00    6.68117e-02
5 43889057 C T  rs13189727      336474  1.25266e+04     4.21457e+02     -1.25947e-02    6.42476e-03     -1.96034e+00    4.99575e-02
5 43889207 A G  rs4516856       336474  6.69113e+05     2.89791e+04     1.96508e-02     1.14359e-02     1.71835e+00     8.57338e-02
5 43889333 G T  rs114787943     336474  3.03587e+03     5.56253e+01     -1.58425e-02    1.28856e-02     -1.22947e+00    2.18896e-01
[  mturchin@node1149  ~]$zcat $HOME/data/mturchin/Neale2017/Vs2/50.assoc.tsv.gz | sed 's/:/ /g' | awk '{ print $5 "\t" $3 "\t" $4 "\t" $9 "\t" $10 "\t" $12 }' | head -n 10
ytx     nCompleteSamples        AC      pval
rs13184706      C       T       -1.33660e-02    6.48300e-03     3.92375e-02
rs58824264      C       T       -1.88438e-02    1.44299e-02     1.91592e-01
rs72762387      T       C       -9.45691e-03    5.62698e-03     9.28345e-02
rs115032754     C       T       1.06178e-03     6.29484e-03     8.66053e-01
rs147555725     C       G       -4.19957e-03    2.06522e-02     8.38864e-01
rs13185925      G       C       -8.69909e-03    4.74596e-03     6.68117e-02
rs13189727      C       T       -1.25947e-02    6.42476e-03     4.99575e-02
rs4516856       A       G       1.96508e-02     1.14359e-02     8.57338e-02
rs114787943     G       T       -1.58425e-02    1.28856e-02     2.18896e-01
[  mturchin@node1149  ~]$zcat $HOME/data/mturchin/Neale2017/Vs2/50.assoc.tsv.gz | awk '{ print $6 "\t" $7 "\t" $9 }' | head -n 10
beta    se      pval
-1.33660e-02    6.48300e-03     3.92375e-02
-1.88438e-02    1.44299e-02     1.91592e-01
-9.45691e-03    5.62698e-03     9.28345e-02
1.06178e-03     6.29484e-03     8.66053e-01
-4.19957e-03    2.06522e-02     8.38864e-01
-8.69909e-03    4.74596e-03     6.68117e-02
-1.25947e-02    6.42476e-03     4.99575e-02
1.96508e-02     1.14359e-02     8.57338e-02
-1.58425e-02    1.28856e-02     2.18896e-01
[  mturchin@node1308  ~]$zcat $HOME/data/mturchin/Loh2017/body_HEIGHTz.sumstats.gz | awk '{ if ($4 != $6) { $8 = -1 * $8; val1 = $4; val2 = $5; $4 = $val2; $5 = $val1; }; print $1 "\t" $4 "\t" $5 "\t" $8 "\t" $9 "\t" $10 }' | head -n 10
SNP     SNP CHR POS A1 A2 REF EAF 0 se P N INFO SNP CHR POS SNP CHR POS A1 A2 REF EAF 0 se P N INFO A2 REF EAF 0 se P N INFO    0       se      P
rs10399793      T       C       -0.00104841     0.00272961      8.6E-01
rs2462492       C       T       -0.00505746     0.00270391      1.1E-01
rs3107975       T       C       0.00918983      0.015069        6.3E-01
rs74447903      T       C       -0.00814975     0.0335016       7.8E-01
1:70728_C_T     C       T       -0.0222717      0.0272397       5.9E-01
rs2462495       A       G       -0.0429876      0.0346811       3.6E-01
rs114608975     T       C       -0.000955315    0.00432242      5.1E-01
rs6702460       G       T       -0.00371767     0.00266258      3.6E-01
rs8179466       C       T       0.00195976      0.00525025      7.0E-01
[  mturchin@node1308  ~]$zcat $HOME/data/mturchin/Loh2017/body_HEIGHTz.sumstats.gz | awk '{ if ($4 != $6) { $8 = -1 * $8; val1 = $4; val2 = $5; $4 = val2; $5 = val1; }; print $1 "\t" $4 "\t" $5 "\t" $8 "\t" $9 "\t" $10 }' | head -n 10
SNP     A2      A1      0       se      P
rs10399793      T       C       -0.00104841     0.00272961      8.6E-01
rs2462492       C       T       -0.00505746     0.00270391      1.1E-01
rs3107975       T       C       0.00918983      0.015069        6.3E-01
rs74447903      T       C       -0.00814975     0.0335016       7.8E-01
1:70728_C_T     C       T       -0.0222717      0.0272397       5.9E-01
rs2462495       A       G       -0.0429876      0.0346811       3.6E-01
rs114608975     T       C       -0.000955315    0.00432242      5.1E-01
rs6702460       G       T       -0.00371767     0.00266258      3.6E-01
rs8179466       C       T       0.00195976      0.00525025      7.0E-01
[  mturchin@node1308  ~]$zcat $HOME/data/mturchin/Loh2017/body_HEIGHTz.sumstats.gz | awk '{ if ($4 != $6) { $8 = -1 * $8; print $4 "\t" $5; val1 = $4; val2 = $5; $4 = val2; $5 = val1; print $4 "\t" $5; }; print $1 "\t" $4 "\t" $5 "\t" $8 "\t" $9 "\t" $10 }' | head -n 10
A1      A2
A2      A1
SNP     A2      A1      0       se      P
rs10399793      T       C       -0.00104841     0.00272961      8.6E-01
rs2462492       C       T       -0.00505746     0.00270391      1.1E-01
rs3107975       T       C       0.00918983      0.015069        6.3E-01
rs74447903      T       C       -0.00814975     0.0335016       7.8E-01
1:70728_C_T     C       T       -0.0222717      0.0272397       5.9E-01
rs2462495       A       G       -0.0429876      0.0346811       3.6E-01
rs114608975     T       C       -0.000955315    0.00432242      5.1E-01
[  mturchin@node1308  ~]$zcat /users/mturchin/data/1000G/mturchin20/Analyses/ShiftAE/Vs1/SubFiles/PerBetaFiles/CEUGBRTSIESNYRIESN.chr22.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.SNPs.noEURfix.edit.wMeanInfo.Loh2017.Height.Inc.frq.gz | head -n 10
rs1000112 G C 0.899 0.8956 0.93925 0.93458 0.9171075 -0.0181074999999999 -0.0215075 0.0221425000000001 0.0174725 G 2
rs1000205 A T 0.4747 0.478 0.5327 0.472 0.48935 -0.0146499999999999 -0.01135 0.04335 -0.01735 A 1
rs1000427 A G 0.1263 0.08242 0.1355 0.06542 0.10241 0.02389 -0.01999 0.03309 -0.03699 A 1
rs1000470 C A 0.899 0.8516 0.92991 0.9486 0.9072775 -0.00827749999999994 -0.0556774999999999 0.0226325000000001 0.0413225 C 2
rs1000539 G A 0.7424 0.6813 0.7103 0.6682 0.70055 0.0418499999999999 -0.01925 0.00975000000000004 -0.03235 G 2
rs1000815 A G 0.4798 0.4725 0.4112 0.4112 0.443675 0.036125 0.028825 -0.032475 -0.032475 A 1
rs10009 G A 0.404 0.4341 0.4346 0.4159 0.42215 -0.01815 0.01195 0.01245 -0.00625000000000003 G 1
rs1001008 G A 0.5 0.489 0.472 0.5421 0.500775 -0.00077499999999997 -0.011775 -0.028775 0.0413250000000001 G 1
rs1001009 C T 0.5 0.4835 0.472 0.5421 0.4994 0.000599999999999989 -0.0159 -0.0274 0.0427 C 1
rs1001021 G A 0.9596 0.94505 0.93925 0.9486 0.948125 0.011475 -0.00307500000000005 -0.00887499999999997 0.000475000000000003 G 2
[  mturchin@node1308  ~]$zcat /users/mturchin/data/1000G/mturchin20/Analyses/ShiftAE/Vs1/SubFiles/PerBetaFiles/CEUGBRTSIESNYRIESN.chr22.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.SNPs.noEURfix.edit.wMeanInfo.Loh2017.Height.Inc.frq.gz | wc        
 105873 1482222 11583858
 [  mturchin@node1308  ~]$zcat /users/mturchin/data/1000G/mturchin20/Analyses/ShiftAE/Vs1/CEUGBRTSIESNYRIESN.chr${chr1}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.SNPs.noEURfix.edit.wMeanInfo.frq.gz | grep rs1000112
 rs1000112 C G 0.101 0.1044 0.06075 0.06542 0.0828925 0.0181075 0.0215075 0.0221425 0.0174725
 [  mturchin@node1308  ~]$zcat /users/mturchin/data/1000G/mturchin20/Analyses/ShiftAE/Vs1/CEUGBRTSIESNYRIESN.chr${chr1}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.SNPs.noEURfix.edit.wMeanInfo.frq.gz | grep rs1000470
 rs1000470 A C 0.101 0.1484 0.07009 0.0514 0.0927225 0.00827750000000001 0.0556775 0.0226325 0.0413225
 [  mturchin@node1308  ~]$zcat /users/mturchin/data/1000G/mturchin20/Analyses/ShiftAE/Vs1/CEUGBRTSIESNYRIESN.chr${chr1}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.SNPs.noEURfix.edit.wMeanInfo.frq.gz | grep rs1001021
 rs1001021 A G 0.0404 0.05495 0.06075 0.0514 0.051875 0.011475 0.003075 0.008875 0.000474999999999996
 [  mturchin@node1308  ~]$zcat /users/mturchin/data/1000G/mturchin20/Analyses/ShiftAE/Vs1/SubFiles/PerBetaFiles/CEUGBRTSIESNYRIESN.chr22.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.SNPs.noEURfix.edit.wMeanInfo.Loh2017.Height.Inc.frq.gz | tail -n 10
 rs998937 A G 0.4091 0.3901 0.3411 0.3925 0.3832 0.0259 0.00690000000000002 -0.0421 0.00930000000000003 A 1
 rs999223 G A 0.2172 0.2582 0.3458 0.2383 0.264875 -0.047675 -0.00667499999999999 0.080925 -0.026575 G 1
 rs999224 G A 0.2222 0.2527 0.3458 0.2383 0.26475 -0.04255 -0.01205 0.08105 -0.02645 G 1
 rs999379 T G 0.4848 0.478 0.5327 0.4533 0.4872 -0.00239999999999996 -0.00919999999999999 0.0455 -0.0339 T 1
 rs999458 C A 0.6818 0.7253 0.7243 0.7243 0.713925 -0.0321250000000001 0.0113749999999999 0.010375 0.010375 C 2
 rs999540 A G 0.1717 0.1374 0.1682 0.1682 0.161375 0.010325 -0.023975 0.006825 0.006825 A 1
 rs999699 T C 0.0404 0.1044 0.04206 0.07944 0.066575 -0.026175 0.037825 -0.024515 0.012865 T 1
 rs9997 C T 0.3333 0.3132 0.3458 0.3224 0.328675 0.00462499999999999 -0.015475 0.017125 -0.00627499999999998 C 1
 rs999719 T A 0.7525 0.6484 0.7243 0.6869 0.703025 0.0494749999999999 -0.054625 0.021275 -0.0161250000000001 T 2
 rs999927 C A 0.6212 0.6538 0.5888 0.6776 0.63535 -0.01415 0.0184500000000001 -0.04655 0.04225 C 2
 [  mturchin@node1308  ~]$zcat /users/mturchin/data/1000G/mturchin20/Analyses/ShiftAE/Vs1/CEUGBRTSIESNYRIESN.chr${chr1}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.SNPs.noEURfix.edit.wMeanInfo.frq.gz | grep rs999458
 rs999458 A C 0.3182 0.2747 0.2757 0.2757 0.286075 0.032125 0.011375 0.010375 0.010375
 [  mturchin@node1308  ~]$zcat /users/mturchin/data/1000G/mturchin20/Analyses/ShiftAE/Vs1/CEUGBRTSIESNYRIESN.chr${chr1}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.SNPs.noEURfix.edit.wMeanInfo.frq.gz | grep rs999719
 rs999719 A T 0.2475 0.3516 0.2757 0.3131 0.296975 0.049475 0.054625 0.021275 0.016125
 [  mturchin@node1308  ~]$zcat /users/mturchin/data/1000G/mturchin20/Analyses/ShiftAE/Vs1/CEUGBRTSIESNYRIESN.chr${chr1}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.SNPs.noEURfix.edit.wMeanInfo.frq.gz | grep rs999927
 rs999927 A C 0.3788 0.3462 0.4112 0.3224 0.36465 0.01415 0.01845 0.04655 0.04225
[  mturchin@node1321  ~]$zcat /users/mturchin/data/1000G/mturchin20/Analyses/ShiftAE/Vs1/CEUGBRTSIESNYRIESN.chr*.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.SNPs.noEURfix.edit.wMeanInfo.frq.gz | awk '{ print $1 }' | sort | wc
8624540 8624540 92623539
[  mturchin@node1321  ~]$zcat /users/mturchin/data/1000G/mturchin20/Analyses/ShiftAE/Vs1/CEUGBRTSIESNYRIESN.chr*.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.SNPs.noEURfix.edit.wMeanInfo.frq.gz | awk '{ print $1 }' | sort | uniq | wc
8623143 8623143 92607576
[  mturchin@node1321  ~]$zcat /users/mturchin/data/1000G/mturchin20/Analyses/ShiftAE/Vs1/CEUGBRTSIESNYRIESN.chr*.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.SNPs.noEURfix.edit.wMeanInfo.frq.gz | awk '{ print $1 }' | sort | uniq -u | wc
8623090 8623090 92606975
(MultiEthnicGWAS) [  mturchin@login003  ~/LabMisc/RamachandranLab/ShiftAE]$zcat /users/mturchin/data/1000G/mturchin20/Analyses/ShiftAE/Vs1/CEUGBRTSIESNYRIESN.chr${chr1}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.SNPs.frq.gz | awk '{ print $1 }' | sort | uniq -c | sort -r -g -k 1,1 | head -n 10
     64 rs9442277
     64 rs78127415
     64 rs587708826
     64 rs587648100
     64 rs551178717
     64 rs549010975
     64 rs548431711
     64 rs533758591
     64 rs202107783
     64 rs201680403
[  mturchin@node1321  ~]$zcat /users/mturchin/data/1000G/mturchin20/Analyses/ShiftAE/Vs1/CEUGBRTSIESNYRIESN.chr${chr1}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.SNPs.frq.gz | grep rs533758591 | head -n 10
rs533758591     0       T       0       0       T       0       0       T       0       0       T       0       0       T       0       0       T       0
rs533758591     0       T       0       0       T       0       0       T       0       0       T       0       0       T       0       0       T       0
rs533758591     0       T       0       0       T       0       0       T       0       0       T       0       A       T       0.00463 0       T       0
rs533758591     0       T       0       0       T       0       0       T       0       0       T       0       A       T       0.00463 0       T       0
rs533758591     0       T       0       0       T       0       0       T       0       A       T       0.004673        0       T       0       0       T       0
rs533758591     0       T       0       0       T       0       0       T       0       A       T       0.004673        0       T       0       0       T       0
rs533758591     0       T       0       0       T       0       0       T       0       A       T       0.004673        A       T       0.00463 0       T       0
rs533758591     0       T       0       0       T       0       0       T       0       A       T       0.004673        A       T       0.00463 0       T       0
rs533758591     0       T       0       0       T       0       0       T       0       0       T       0       0       T       0       0       T       0
rs533758591     0       T       0       0       T       0       0       T       0       0       T       0       0       T       0       0       T       0
(MultiEthnicGWAS) [  mturchin@login003  ~/LabMisc/RamachandranLab/ShiftAE]$zcat /users/mturchin/data/1000G/mturchin20/Analyses/ShiftAE/Vs1/CEUGBRTSIESNYRIESN.chr${chr1}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.SNPs.frq.gz | awk '{ print $1 }' | sort | uniq -c | sort -r -g -k 1,1 | head -n 10
      1 rs999992   
      1 rs999974
      1 rs999933
      1 rs999932
      1 rs999923
      1 rs999918
      1 rs999917
      1 rs999915
      1 rs999897
      1 rs999894
[  mturchin@node1308  ~]$zcat /users/mturchin/data/1000G/mturchin20/Analyses/ShiftAE/Vs1/SubFiles/PerAFFiles/*frq.gz | head -n 10
rs1001219 T C 0.0101 0.02198 0.004673 0.009346 0.01152475 0.00142475 0.01045525 0.00685175 0.00217875
rs1001246 C T 0.0101 0.005495 0.009346 0.01402 0.00974025 0.000359750000000001 0.00424525 0.000394249999999999 0.00427975
rs1002194 T C 0.0303 0.01099 0.004673 0.004673 0.012659 0.017641 0.001669 0.007986 0.007986
rs1002383 G T 0.005051 0.005495 0.004673 0.009346 0.00614125 0.00109025 0.00064625 0.00146825 0.00320475
rs1005852 C T 0.0202 0.01099 0.01402 0.004673 0.01247075 0.00772925 0.00148075 0.00154925 0.00779775
rs1007860 T G 0.005051 0.01648 0.009346 0.009346 0.01005575 0.00500475 0.00642425 0.00070975 0.00070975
rs10082122 C T 0.01515 0.01099 0.004673 0.01402 0.01120825 0.00394175 0.00021825 0.00653525 0.00281175
rs10082157 T C 0.01515 0.005495 0.004673 0.009346 0.008666 0.006484 0.003171 0.003993 0.00068
rs1009127 T G 0.0101 0.005495 0.004673 0.004673 0.00623525 0.00386475 0.00074025 0.00156225 0.00156225
rs1010283 A T 0.01515 0.01099 0.009346 0.004673 0.01003975 0.00511025 0.00095025 0.00069375 0.00536675
[  mturchin@node1308  ~]$zcat /users/mturchin/data/1000G/mturchin20/Analyses/ShiftAE/Vs1/SubFiles/PerAFFiles/*frq.gz | awk '{ print $1 }' | sort | wc
8623077 8623077 92606827
[  mturchin@node1308  ~]$zcat /users/mturchin/data/1000G/mturchin20/Analyses/ShiftAE/Vs1/SubFiles/PerAFFiles/*frq.gz | awk '{ print $1 }' | sort | uniq | wc
8623077 8623077 92606827



~~~



