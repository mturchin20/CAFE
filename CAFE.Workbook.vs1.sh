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


#20190321
#Vs1

#Idea -- SNPs undergoing drift should randomly fluctuate to both higher and lower frequencies 'equally' (modulo population demographic misc) from their ancestral state over time. SNPs undergoing weak selective pressures in a given direction should experience less overall variance in their change in frequency since their ancestral state (alleles that are favored get pushed up; if you don't know orientation properly, can just use absolute difference since amount pushed up or amount pushed down are equivalent; looking at overall standard difference amongst 'AF change since ancestral state'; AF affects how much has changed since ancestral, but will be doing matched permutations on ancestral AF so within the final score, that's fine; the idea is that all alleles have similar changes in AF since ancestral state when they're being weakly selected, so they are more similar than random chance alone given a normal distribution variance)

mkdir /users/mturchin/LabMisc/RamachandranLab/InterPath/CAFE
mkdir /users/mturchin/LabMisc/RamachandranLab/InterPath/CAFE/Vs1
mkdir /users/mturchin/LabMisc/RamachandranLab/InterPath/CAFE/Vs1/Analyses

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
#From midway
#scp -p /project/mstephens/mturchin20/gluster/data/external_public_supp/GIANT2014_5/GIANT_HEIGHT_Wood_et_al_2014_publicrelease_HapMapCeuFreq.wUCSCGenomeBrowser_dbSNP130.vs1.txt.gz mturchin@ssh.ccv.brown.edu:/users/mturchin/Data2/GIANT/.

zcat $HOME/data/mturchin/Neale2017/Vs2/50.assoc.tsv.gz | awk '{ if ($9 < 5e-9) { print $0 } } ' | gzip > $HOME/data/mturchin/Neale2017/Vs2/50.assoc.lt5eNeg9.tsv.gz
zcat $HOME/data/mturchin/Neale2017/Vs2/50.assoc.tsv.gz | awk '{ if ($9 < 1e-4) { print $0 } } ' | gzip > $HOME/data/mturchin/Neale2017/Vs2/50.assoc.lt1eNeg4.tsv.gz
zcat $HOME/data/mturchin/Loh2017/body_HEIGHTz.sumstats.gz | awk '{ if ($10 < 5e-9) { print $0 } } ' | gzip > $HOME/data/mturchin/Loh2017/body_HEIGHTz.sumstats.lt5eNeg9.gz
zcat $HOME/data/mturchin/Loh2017/body_HEIGHTz.sumstats.gz | awk '{ if ($10 < 1e-4) { print $0 } } ' | gzip > $HOME/data/mturchin/Loh2017/body_HEIGHTz.sumstats.lt1eNeg4.gz
zcat /users/mturchin/Data2/GIANT/GIANT_HEIGHT_Wood_et_al_2014_publicrelease_HapMapCeuFreq.wUCSCGenomeBrowser_dbSNP130.vs1.txt.gz | awk '{ if ($7 < 5e-9) { print $0 } } ' | gzip > /users/mturchin/Data2/GIANT/GIANT_HEIGHT_Wood_et_al_2014_publicrelease_HapMapCeuFreq.wUCSCGenomeBrowser_dbSNP130.vs1.lt5eNeg9.txt.gz
zcat /users/mturchin/Data2/GIANT/GIANT_HEIGHT_Wood_et_al_2014_publicrelease_HapMapCeuFreq.wUCSCGenomeBrowser_dbSNP130.vs1.txt.gz | awk '{ if ($7 < 5e-8) { print $0 } } ' | gzip > /users/mturchin/Data2/GIANT/GIANT_HEIGHT_Wood_et_al_2014_publicrelease_HapMapCeuFreq.wUCSCGenomeBrowser_dbSNP130.vs1.lt5eNeg8.txt.gz
zcat /users/mturchin/Data2/GIANT/GIANT_HEIGHT_Wood_et_al_2014_publicrelease_HapMapCeuFreq.wUCSCGenomeBrowser_dbSNP130.vs1.txt.gz | awk '{ if ($7 < 1e-4) { print $0 } } ' | gzip > /users/mturchin/Data2/GIANT/GIANT_HEIGHT_Wood_et_al_2014_publicrelease_HapMapCeuFreq.wUCSCGenomeBrowser_dbSNP130.vs1.lt1eNeg4.txt.gz

##zcat $HOME/data/mturchin/Loh2017/body_HEIGHTz.sumstats.lt5eNeg9.gz | R -q -e "library(\"devtools\"); devtools::load_all(\"/users/mturchin/LabMisc/StephensLab/bmass\"); Data1 <- read.table(file('stdin'), header=F); for (i in 1:22) { print(i); Data2 <- Data1[Data1[,2] == i,]; Data3 <- Data2[order(-log10(Data2[,10]), decreasing=TRUE),]; Data3 <- Data3[indephits(-log10(Data2[,10]), Data2[,2], Data2[,3]) == 1,]; print(dim(Data3)); print(head(Data3)); };"
rm $HOME/data/mturchin/Loh2017/body_HEIGHTz.sumstats.lt5eNeg9.bmass.GrdyClmp.txt.gz; zcat $HOME/data/mturchin/Loh2017/body_HEIGHTz.sumstats.lt5eNeg9.gz | R -q -e "library(\"devtools\"); devtools::load_all(\"/users/mturchin/LabMisc/StephensLab/bmass\"); Data1 <- read.table(file('stdin'), header=F); for (i in 1:22) { print(i); Data2 <- Data1[Data1[,2] == i,]; Data3 <- Data2[order(-log10(Data2[,10]), decreasing=TRUE),]; Data3 <- Data3[indephits(-log10(Data2[,10]), Data2[,2], Data2[,3]) == 1,]; write.table(Data3, file=\"$HOME/data/mturchin/Loh2017/body_HEIGHTz.sumstats.lt5eNeg9.bmass.GrdyClmp.txt\", append=TRUE, quote=FALSE, col.names=FALSE, row.names=FALSE); };"; gzip -f $HOME/data/mturchin/Loh2017/body_HEIGHTz.sumstats.lt5eNeg9.bmass.GrdyClmp.txt
rm -f $HOME/data/mturchin/Loh2017/body_HEIGHTz.sumstats.lt1eNeg4.bmass.GrdyClmp.txt.gz; zcat $HOME/data/mturchin/Loh2017/body_HEIGHTz.sumstats.lt1eNeg4.gz | R -q -e "library(\"devtools\"); devtools::load_all(\"/users/mturchin/LabMisc/StephensLab/bmass\"); Data1 <- read.table(file('stdin'), header=F); for (i in 1:22) { print(i); Data2 <- Data1[Data1[,2] == i,]; Data3 <- Data2[order(-log10(Data2[,10]), decreasing=TRUE),]; Data3 <- Data3[indephits(-log10(Data2[,10]), Data2[,2], Data2[,3]) == 1,]; write.table(Data3, file=\"$HOME/data/mturchin/Loh2017/body_HEIGHTz.sumstats.lt1eNeg4.bmass.GrdyClmp.txt\", append=TRUE, quote=FALSE, col.names=FALSE, row.names=FALSE); };"; gzip -f $HOME/data/mturchin/Loh2017/body_HEIGHTz.sumstats.lt1eNeg4.bmass.GrdyClmp.txt
rm -f $HOME/data/mturchin/Neale2017/Vs2/50.assoc.lt5eNeg9.bmass.GrdyClmp.tsv.gz; zcat $HOME/data/mturchin/Neale2017/Vs2/50.assoc.lt5eNeg9.tsv.gz | sed 's/:/ /' | sed 's/:/ /' | R -q -e "library(\"devtools\"); devtools::load_all(\"/users/mturchin/LabMisc/StephensLab/bmass\"); Data1 <- read.table(file('stdin'), header=F); for (i in 1:22) { print(i); Data2 <- Data1[Data1[,1] == i,]; Data3 <- Data2[order(-log10(Data2[,11]), decreasing=TRUE),]; Data3 <- Data3[indephits(-log10(Data2[,11]), Data2[,1], Data2[,2]) == 1,]; write.table(Data3, file=\"$HOME/data/mturchin/Neale2017/Vs2/50.assoc.lt5eNeg9.bmass.GrdyClmp.tsv\", append=TRUE, quote=FALSE, col.names=FALSE, row.names=FALSE); };"; gzip -f $HOME/data/mturchin/Neale2017/Vs2/50.assoc.lt5eNeg9.bmass.GrdyClmp.tsv
rm -f $HOME/data/mturchin/Neale2017/Vs2/50.assoc.lt1eNeg4.bmass.GrdyClmp.tsv.gz; zcat $HOME/data/mturchin/Neale2017/Vs2/50.assoc.lt1eNeg4.tsv.gz | sed 's/:/ /' | sed 's/:/ /' | R -q -e "library(\"devtools\"); devtools::load_all(\"/users/mturchin/LabMisc/StephensLab/bmass\"); Data1 <- read.table(file('stdin'), header=F); for (i in 1:22) { print(i); Data2 <- Data1[Data1[,1] == i,]; Data3 <- Data2[order(-log10(Data2[,11]), decreasing=TRUE),]; Data3 <- Data3[indephits(-log10(Data2[,11]), Data2[,1], Data2[,2]) == 1,]; write.table(Data3, file=\"$HOME/data/mturchin/Neale2017/Vs2/50.assoc.lt1eNeg4.bmass.GrdyClmp.tsv\", append=TRUE, quote=FALSE, col.names=FALSE, row.names=FALSE); };"; gzip -f $HOME/data/mturchin/Neale2017/Vs2/50.assoc.lt1eNeg4.bmass.GrdyClmp.tsv
rm -f /users/mturchin/Data2/GIANT/GIANT_HEIGHT_Wood_et_al_2014_publicrelease_HapMapCeuFreq.wUCSCGenomeBrowser_dbSNP130.vs1.lt5eNeg9.GrdyClmp.tx*; zcat /users/mturchin/Data2/GIANT/GIANT_HEIGHT_Wood_et_al_2014_publicrelease_HapMapCeuFreq.wUCSCGenomeBrowser_dbSNP130.vs1.lt5eNeg9.txt.gz | grep -v NA | grep -v , | sed 's/_/ /g' | R -q -e "library(\"devtools\"); devtools::load_all(\"/users/mturchin/LabMisc/StephensLab/bmass\"); Data1 <- read.table(file('stdin'), header=F); for (i in 1:22) { print(i); Data2 <- Data1[Data1[,9] == i,]; Data3 <- Data2[order(-log10(Data2[,7]), decreasing=TRUE),]; Data3 <- Data3[indephits(-log10(Data2[,7]), Data2[,9], Data2[,10]) == 1,]; write.table(Data3, file=\"/users/mturchin/Data2/GIANT/GIANT_HEIGHT_Wood_et_al_2014_publicrelease_HapMapCeuFreq.wUCSCGenomeBrowser_dbSNP130.vs1.lt5eNeg9.GrdyClmp.txt\", append=TRUE, quote=FALSE, col.names=FALSE, row.names=FALSE); };"; gzip -f /users/mturchin/Data2/GIANT/GIANT_HEIGHT_Wood_et_al_2014_publicrelease_HapMapCeuFreq.wUCSCGenomeBrowser_dbSNP130.vs1.lt5eNeg9.GrdyClmp.txt
rm -f /users/mturchin/Data2/GIANT/GIANT_HEIGHT_Wood_et_al_2014_publicrelease_HapMapCeuFreq.wUCSCGenomeBrowser_dbSNP130.vs1.lt5eNeg8.GrdyClmp.tx*; zcat /users/mturchin/Data2/GIANT/GIANT_HEIGHT_Wood_et_al_2014_publicrelease_HapMapCeuFreq.wUCSCGenomeBrowser_dbSNP130.vs1.lt5eNeg8.txt.gz | grep -v NA | grep -v , | sed 's/_/ /g' | R -q -e "library(\"devtools\"); devtools::load_all(\"/users/mturchin/LabMisc/StephensLab/bmass\"); Data1 <- read.table(file('stdin'), header=F); for (i in 1:22) { print(i); Data2 <- Data1[Data1[,9] == i,]; Data3 <- Data2[order(-log10(Data2[,7]), decreasing=TRUE),]; Data3 <- Data3[indephits(-log10(Data2[,7]), Data2[,9], Data2[,10]) == 1,]; write.table(Data3, file=\"/users/mturchin/Data2/GIANT/GIANT_HEIGHT_Wood_et_al_2014_publicrelease_HapMapCeuFreq.wUCSCGenomeBrowser_dbSNP130.vs1.lt5eNeg8.GrdyClmp.txt\", append=TRUE, quote=FALSE, col.names=FALSE, row.names=FALSE); };"; gzip -f /users/mturchin/Data2/GIANT/GIANT_HEIGHT_Wood_et_al_2014_publicrelease_HapMapCeuFreq.wUCSCGenomeBrowser_dbSNP130.vs1.lt5eNeg8.GrdyClmp.txt
rm -f /users/mturchin/Data2/GIANT/GIANT_HEIGHT_Wood_et_al_2014_publicrelease_HapMapCeuFreq.wUCSCGenomeBrowser_dbSNP130.vs1.lt1eNeg4.GrdyClmp.tx*; zcat /users/mturchin/Data2/GIANT/GIANT_HEIGHT_Wood_et_al_2014_publicrelease_HapMapCeuFreq.wUCSCGenomeBrowser_dbSNP130.vs1.lt1eNeg4.txt.gz | grep -v NA | grep -v , | sed 's/_/ /g' | R -q -e "library(\"devtools\"); devtools::load_all(\"/users/mturchin/LabMisc/StephensLab/bmass\"); Data1 <- read.table(file('stdin'), header=F); for (i in 1:22) { print(i); Data2 <- Data1[Data1[,9] == i,]; Data3 <- Data2[order(-log10(Data2[,7]), decreasing=TRUE),]; Data3 <- Data3[indephits(-log10(Data2[,7]), Data2[,9], Data2[,10]) == 1,]; write.table(Data3, file=\"/users/mturchin/Data2/GIANT/GIANT_HEIGHT_Wood_et_al_2014_publicrelease_HapMapCeuFreq.wUCSCGenomeBrowser_dbSNP130.vs1.lt1eNeg4.GrdyClmp.txt\", append=TRUE, quote=FALSE, col.names=FALSE, row.names=FALSE); };"; gzip -f /users/mturchin/Data2/GIANT/GIANT_HEIGHT_Wood_et_al_2014_publicrelease_HapMapCeuFreq.wUCSCGenomeBrowser_dbSNP130.vs1.lt1eNeg4.GrdyClmp.txt

zcat /users/mturchin/data/mturchin/Loh2017/body_HEIGHTz.sumstats.lt5eNeg9.bmass.GrdyClmp.txt.gz | awk '{ print $1 }' | gzip > /users/mturchin/data/mturchin/Loh2017/body_HEIGHTz.sumstats.lt5eNeg9.bmass.GrdyClmp.rsIDs.gz 
zcat /users/mturchin/data/mturchin/Loh2017/body_HEIGHTz.sumstats.lt1eNeg4.bmass.GrdyClmp.txt.gz | awk '{ print $1 }' | gzip > /users/mturchin/data/mturchin/Loh2017/body_HEIGHTz.sumstats.lt1eNeg4.bmass.GrdyClmp.rsIDs.gz 
zcat /users/mturchin/data/mturchin//Neale2017/Vs2/50.assoc.lt5eNeg9.bmass.GrdyClmp.tsv.gz | awk '{ print $4 }' | gzip > /users/mturchin/data/mturchin//Neale2017/Vs2/50.assoc.lt5eNeg9.bmass.GrdyClmp.rsIDs.gz
zcat /users/mturchin/data/mturchin//Neale2017/Vs2/50.assoc.lt1eNeg4.bmass.GrdyClmp.tsv.gz | awk '{ print $4 }' | gzip > /users/mturchin/data/mturchin//Neale2017/Vs2/50.assoc.lt1eNeg4.bmass.GrdyClmp.rsIDs.gz
zcat /users/mturchin/Data2/GIANT/GIANT_HEIGHT_Wood_et_al_2014_publicrelease_HapMapCeuFreq.wUCSCGenomeBrowser_dbSNP130.vs1.lt5eNeg9.GrdyClmp.txt.gz | awk '{ print $1 }' | gzip > /users/mturchin/Data2/GIANT/GIANT_HEIGHT_Wood_et_al_2014_publicrelease_HapMapCeuFreq.wUCSCGenomeBrowser_dbSNP130.vs1.lt5eNeg9.GrdyClmp.rsIDs.gz 
zcat /users/mturchin/Data2/GIANT/GIANT_HEIGHT_Wood_et_al_2014_publicrelease_HapMapCeuFreq.wUCSCGenomeBrowser_dbSNP130.vs1.lt5eNeg8.GrdyClmp.txt.gz | awk '{ print $1 }' | gzip > /users/mturchin/Data2/GIANT/GIANT_HEIGHT_Wood_et_al_2014_publicrelease_HapMapCeuFreq.wUCSCGenomeBrowser_dbSNP130.vs1.lt5eNeg8.GrdyClmp.rsIDs.gz 
zcat /users/mturchin/Data2/GIANT/GIANT_HEIGHT_Wood_et_al_2014_publicrelease_HapMapCeuFreq.wUCSCGenomeBrowser_dbSNP130.vs1.lt1eNeg4.GrdyClmp.txt.gz | awk '{ print $1 }' | gzip > /users/mturchin/Data2/GIANT/GIANT_HEIGHT_Wood_et_al_2014_publicrelease_HapMapCeuFreq.wUCSCGenomeBrowser_dbSNP130.vs1.lt1eNeg4.GrdyClmp.rsIDs.gz 

#/users/mturchin/LabMisc/HirschhornLab/SohailRspnd/ukb_chrAll_v2.British.QCed.reqDrop.QCed.dropRltvs.PCAdrop.Height.Trans.ADD.assoc.linear.LohpVals.ppr.clumped.gz
#/users/mturchin/LabMisc/HirschhornLab/SohailRspnd/ukb_chrAll_v2.British.QCed.reqDrop.QCed.dropRltvs.PCAdrop.Height.Trans.ADD.assoc.linear.LohpVals.loose.clumped.gz
#/users/mturchin/LabMisc/HirschhornLab/SohailRspnd/ukb_chrAll_v2.British.QCed.reqDrop.QCed.dropRltvs.PCAdrop.Height.Trans.ADD.assoc.linear.BNpVals.ppr.clumped.gz
#/users/mturchin/LabMisc/HirschhornLab/SohailRspnd/ukb_chrAll_v2.British.QCed.reqDrop.QCed.dropRltvs.PCAdrop.Height.Trans.ADD.assoc.linear.BNpVals.loose.clumped.gz
##/users/mturchin/LabMisc/HirschhornLab/SohailRspnd/ukb_chrAll_v2.British.QCed.reqDrop.QCed.dropRltvs.PCAdrop.Height.Trans.ADD.assoc.linear.GIANTpVals.ppr.clumped.gz
##/users/mturchin/LabMisc/HirschhornLab/SohailRspnd/ukb_chrAll_v2.British.QCed.reqDrop.QCed.dropRltvs.PCAdrop.Height.Trans.ADD.assoc.linear.GIANTpVals.loose.clumped.gz
#/users/mturchin/LabMisc/HirschhornLab/20180601_Wood2014_HeightGWASSNPs.txt

mkdir /users/mturchin/data/1000G/mturchin20/Analyses
mkdir /users/mturchin/data/1000G/mturchin20/Analyses/CAFE
mkdir /users/mturchin/data/1000G/mturchin20/Analyses/CAFE/Vs1

#NOTE -- all 1000G files are on the same strand, and since pop subsets were extracted from the main source AllPops files, they should all start with the same set of SNPs at least (read: http://www.internationalgenome.org/faq/what-strand-are-variants-your-vcf-file/)

for chr1 in `echo {1..22}`; do
	echo $chr1

#	join <(zcat /users/mturchin/data/1000G/subsets/CEU/AFs/CEU.chr${chr1}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.SNPs.frq.gz | awk '{ if ($6 > 125) { print $0 } } ' | awk '{ print $2 "\t" $3 "\t" $4 "\t" $5}' | sort -k 1,1) <(zcat /users/mturchin/data/1000G/subsets/GBR/AFs/GBR.chr${chr1}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.SNPs.frq.gz | awk '{ if ($6 > 125) { print $0 } } ' | awk '{ print $2 "\t" $3 "\t" $4 "\t" $5}' | sort -k 1,1) | perl -lane 'if ($F[1] eq "0") { if ($F[4] eq "0") { if ($F[2] eq $F[5]) { print join("\t", @F); } } else { if ($F[2] eq $F[5]) { print join("\t", @F); } elsif ($F[2] eq $F[4]) { ($F[4], $F[5]) = ($F[5], $F[4]); $F[6] = 1 - $F[6]; print join("\t", @F); } else { $PH1 = 1; } } } else { if ($F[4] eq "0") { if ($F[1] eq $F[5]) { ($F[4], $F[5]) = ($F[5], $F[4]); $F[6] = 1 - $F[6]; print join("\t", @F); } elsif ($F[2] eq $F[5]) { print join("\t", @F); } else { $PH1 = 1; } } else { if (($F[1] eq $F[4]) && ($F[2] eq $F[5])) { print join("\t", @F); } elsif (($F[1] eq $F[5]) && ($F[2] eq $F[4])) { ($F[4], $F[5]) = ($F[5], $F[4]); $F[6] = 1 - $F[6]; print join("\t", @F); } else { $PH1 = 1; } } }' | \ 
#	join - <(zcat /users/mturchin/data/1000G/subsets/TSI/AFs/TSI.chr${chr1}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.SNPs.frq.gz | awk '{ if ($6 > 125) { print $0 } } ' | awk '{ print $2 "\t" $3 "\t" $4 "\t" $5}' | sort -k 1,1) | perl -lane 'if ($F[1] eq "0") { if ($F[7] eq "0") { if ($F[2] eq $F[8]) { print join("\t", @F); } } else { if ($F[2] eq $F[8]) { print join("\t", @F); } elsif ($F[2] eq $F[7]) { ($F[7], $F[8]) = ($F[8], $F[7]); $F[9] = 1 - $F[9]; print join("\t", @F); } else { $PH1 = 1; } } } else { if ($F[7] eq "0") { if ($F[1] eq $F[8]) { ($F[7], $F[8]) = ($F[8], $F[7]); $F[9] = 1 - $F[9]; print join("\t", @F); } elsif ($F[2] eq $F[8]) { print join("\t", @F); } else { $PH1 = 1; } } else { if (($F[1] eq $F[7]) && ($F[2] eq $F[8])) { print join("\t", @F); } elsif (($F[1] eq $F[8]) && ($F[2] eq $F[7])) { ($F[7], $F[8]) = ($F[8], $F[7]); $F[9] = 1 - $F[9]; print join("\t", @F); } else { $PH1 = 1; } } }' | \
#	join - <(zcat /users/mturchin/data/1000G/subsets/IBS/AFs/IBS.chr${chr1}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.SNPs.frq.gz | awk '{ if ($6 > 125) { print $0 } } ' | awk '{ print $2 "\t" $3 "\t" $4 "\t" $5}' | sort -k 1,1) | perl -lane 'if ($F[1] eq "0") { if ($F[10] eq "0") { if ($F[2] eq $F[11]) { print join("\t", @F); } } else { if ($F[2] eq $F[11]) { print join("\t", @F); } elsif ($F[2] eq $F[10]) { ($F[10], $F[11]) = ($F[11], $F[10]); $F[12] = 1 - $F[12]; print join("\t", @F); } else { $PH1 = 1; } } } else { if ($F[10] eq "0") { if ($F[1] eq $F[11]) { ($F[10], $F[11]) = ($F[11], $F[10]); $F[12] = 1 - $F[12]; print join("\t", @F); } elsif ($F[2] eq $F[11]) { print join("\t", @F); } else { $PH1 = 1; } } else { if (($F[1] eq $F[10]) && ($F[2] eq $F[11])) { print join("\t", @F); } elsif (($F[1] eq $F[11]) && ($F[2] eq $F[10])) { ($F[10], $F[11]) = ($F[11], $F[10]); $F[12] = 1 - $F[12]; print join("\t", @F); } else { $PH1 = 1; } } }' | \ 
#	join - <(zcat /users/mturchin/data/1000G/subsets/YRI/AFs/YRI.chr${chr1}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.SNPs.frq.gz | awk '{ if ($6 > 125) { print $0 } } ' | awk '{ print $2 "\t" $3 "\t" $4 "\t" $5}' | sort -k 1,1) | perl -lane 'if ($F[1] eq "0") { if ($F[13] eq "0") { if ($F[2] eq $F[14]) { print join("\t", @F); } } else { if ($F[2] eq $F[14]) { print join("\t", @F); } elsif ($F[2] eq $F[13]) { ($F[13], $F[14]) = ($F[14], $F[13]); $F[15] = 1 - $F[15]; print join("\t", @F); } else { $PH1 = 1; } } } else { if ($F[13] eq "0") { if ($F[1] eq $F[14]) { ($F[13], $F[14]) = ($F[14], $F[13]); $F[15] = 1 - $F[15]; print join("\t", @F); } elsif ($F[2] eq $F[14]) { print join("\t", @F); } else { $PH1 = 1; } } else { if (($F[1] eq $F[13]) && ($F[2] eq $F[14])) { print join("\t", @F); } elsif (($F[1] eq $F[14]) && ($F[2] eq $F[13])) { ($F[13], $F[14]) = ($F[14], $F[13]); $F[15] = 1 - $F[15]; print join("\t", @F); } else { $PH1 = 1; } } }' | \ 
#	join - <(zcat /users/mturchin/data/1000G/subsets/ESN/AFs/ESN.chr${chr1}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.SNPs.frq.gz | awk '{ if ($6 > 125) { print $0 } } ' | awk '{ print $2 "\t" $3 "\t" $4 "\t" $5}' | sort -k 1,1) | perl -lane 'if ($F[1] eq "0") { if ($F[16] eq "0") { if ($F[2] eq $F[17]) { print join("\t", @F); } } else { if ($F[2] eq $F[17]) { print join("\t", @F); } elsif ($F[2] eq $F[16]) { ($F[16], $F[17]) = ($F[17], $F[16]); $F[18] = 1 - $F[18]; print join("\t", @F); } else { $PH1 = 1; } } } else { if ($F[16] eq "0") { if ($F[1] eq $F[17]) { ($F[16], $F[17]) = ($F[17], $F[16]); $F[18] = 1 - $F[18]; print join("\t", @F); } elsif ($F[2] eq $F[17]) { print join("\t", @F); } else { $PH1 = 1; } } else { if (($F[1] eq $F[16]) && ($F[2] eq $F[17])) { print join("\t", @F); } elsif (($F[1] eq $F[17]) && ($F[2] eq $F[16])) { ($F[16], $F[17]) = ($F[17], $F[16]); $F[18] = 1 - $F[18]; print join("\t", @F); } else { $PH1 = 1; } } }' | gzip > /users/mturchin/data/1000G/mturchin20/Analyses/CAFE/Vs1/CEUGBRTSIESNYRIESN.chr${chr1}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.SNPs.frq.gz
#	zcat /users/mturchin/data/1000G/mturchin20/Analyses/CAFE/Vs1/CEUGBRTSIESNYRIESN.chr${chr1}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.SNPs.frq.gz | awk '{ if (($4 > 0) && ($7 > 0) && ($10 > 0) && ($13 > 0)) { print $0 } } ' | gzip > /users/mturchin/data/1000G/mturchin20/Analyses/CAFE/Vs1/CEUGBRTSIESNYRIESN.chr${chr1}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.SNPs.noEURfix.frq.gz
	zcat /users/mturchin/data/1000G/mturchin20/Analyses/CAFE/Vs1/CEUGBRTSIESNYRIESN.chr${chr1}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.SNPs.noEURfix.frq.gz | grep -v MAF | R -q -e "Data1 <- read.table(file('stdin'), header=F); EurMean <- apply(Data1[,c(4,7,10,13)], 1, function(x) { return(mean(x)); }); CEUDiff <- abs(Data1[,4] - EurMean); GBRDiff <- abs(Data1[,7] - EurMean); TSIDiff <- abs(Data1[,10] - EurMean); IBSDiff <- abs(Data1[,13] - EurMean); write.table(cbind(Data1[,c(1,2,3,4,7,10,13)], EurMean, CEUDiff, GBRDiff, TSIDiff, IBSDiff), quote=FALSE, col.names=FALSE, row.names=FALSE);" | grep -v \> | gzip > /users/mturchin/data/1000G/mturchin20/Analyses/CAFE/Vs1/CEUGBRTSIESNYRIESN.chr${chr1}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.SNPs.noEURfix.edit.wMeanInfo.frq.gz

done

mkdir /users/mturchin/data/1000G/mturchin20/Analyses/CAFE/Vs1/SubFiles
mkdir /users/mturchin/data/1000G/mturchin20/Analyses/CAFE/Vs1/SubFiles/PerAFFiles

rm -f /users/mturchin/data/1000G/mturchin20/Analyses/CAFE/Vs1/SubFiles/PerAFFiles/*; for chr1 in `echo {1..22}`; do
	echo $chr1

	zcat /users/mturchin/data/1000G/mturchin20/Analyses/CAFE/Vs1/CEUGBRTSIESNYRIESN.chr${chr1}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.SNPs.noEURfix.edit.wMeanInfo.frq.gz | R -q -e "Data1 <- read.table(file('stdin'), header=F); for (i in 1:nrow(Data1)) { AncAF <- round(Data1[i,8], digits=2); write.table(Data1[i,], file=paste(\"/users/mturchin/data/1000G/mturchin20/Analyses/CAFE/Vs1/SubFiles/PerAFFiles/CEUGBRTSIESNYRIESN.chrAll.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.SNPs.noEURfix.edit.wMeanInfo.permPrep.AncAF_\", as.character(AncAF), \".frq\", sep=\"\"), append=TRUE, quote=FALSE, col.names=FALSE, row.names=FALSE);};"

done; gzip -f /users/mturchin/data/1000G/mturchin20/Analyses/CAFE/Vs1/SubFiles/PerAFFiles/*

mkdir /users/mturchin/LabMisc/RamachandranLab/InterPath/CAFE/Data
mkdir /users/mturchin/LabMisc/RamachandranLab/InterPath/CAFE/Data/GWASsnps

ln -s /users/mturchin/data/mturchin/Loh2017/body_HEIGHTz.sumstats.lt5eNeg9.bmass.GrdyClmp.rsIDs.gz /users/mturchin/LabMisc/RamachandranLab/InterPath/CAFE/Data/GWASsnps/Loh2017.Height.lt5eNeg9.bmass.GrdyClmp.rsIDs.gz
ln -s /users/mturchin/data/mturchin/Loh2017/body_HEIGHTz.sumstats.lt1eNeg4.bmass.GrdyClmp.rsIDs.gz /users/mturchin/LabMisc/RamachandranLab/InterPath/CAFE/Data/GWASsnps/Loh2017.Height.lt1eNeg4.bmass.GrdyClmp.rsIDs.gz
ln -s /users/mturchin/data/mturchin//Neale2017/Vs2/50.assoc.lt5eNeg9.bmass.GrdyClmp.rsIDs.gz /users/mturchin/LabMisc/RamachandranLab/InterPath/CAFE/Data/GWASsnps/Neale2017.Height.lt5eNeg9.bmass.GrdyClmp.rsIDs.gz
ln -s /users/mturchin/data/mturchin//Neale2017/Vs2/50.assoc.lt1eNeg4.bmass.GrdyClmp.rsIDs.gz /users/mturchin/LabMisc/RamachandranLab/InterPath/CAFE/Data/GWASsnps/Neale2017.Height.lt1eNeg4.bmass.GrdyClmp.rsIDs.gz
ln -s /users/mturchin/Data2/GIANT/GIANT_HEIGHT_Wood_et_al_2014_publicrelease_HapMapCeuFreq.wUCSCGenomeBrowser_dbSNP130.vs1.lt5eNeg9.GrdyClmp.rsIDs.gz /users/mturchin/LabMisc/RamachandranLab/InterPath/CAFE/Data/GWASsnps/GIANT2014_5.Height.lt5eNeg9.bmass.GrdyClmp.rsIDs.gz
ln -s /users/mturchin/Data2/GIANT/GIANT_HEIGHT_Wood_et_al_2014_publicrelease_HapMapCeuFreq.wUCSCGenomeBrowser_dbSNP130.vs1.lt5eNeg8.GrdyClmp.rsIDs.gz /users/mturchin/LabMisc/RamachandranLab/InterPath/CAFE/Data/GWASsnps/GIANT2014_5.Height.lt5eNeg8.bmass.GrdyClmp.rsIDs.gz
ln -s /users/mturchin/Data2/GIANT/GIANT_HEIGHT_Wood_et_al_2014_publicrelease_HapMapCeuFreq.wUCSCGenomeBrowser_dbSNP130.vs1.lt1eNeg4.GrdyClmp.rsIDs.gz /users/mturchin/LabMisc/RamachandranLab/InterPath/CAFE/Data/GWASsnps/GIANT2014_5.Height.lt1eNeg4.bmass.GrdyClmp.rsIDs.gz
	
for i in `cat <(echo "Loh2017 Neale2017 GIANT2014_5" | perl -lane 'print join("\n", @F);') | head -n 1 | tail -n 1`; do
	for j in `cat <(echo "Height" | perl -lane 'print join("\n", @F);') | head -n 1 | tail -n 1`; do
		for k in `cat <(echo "lt5eNeg9 lt1eNeg4 lt5eNeg8" | perl -lane 'print join("\n", @F);') | head -n 1 | tail -n 1`; do
			echo $i $j $k

			join <(zcat /users/mturchin/LabMisc/RamachandranLab/InterPath/CAFE/Data/GWASsnps/$i.$j.$k.bmass.GrdyClmp.rsIDs.gz | sort) <(zcat /users/mturchin/data/1000G/mturchin20/Analyses/CAFE/Vs1/CEUGBRTSIESNYRIESN.chr*.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.SNPs.noEURfix.edit.wMeanInfo.frq.gz | sort -k 1,1) | gzip > /users/mturchin/LabMisc/RamachandranLab/InterPath/CAFE/Data/GWASsnps/$i.$j.$k.bmass.GrdyClmp.rsIDs.w1000GInfo.txt.gz

		done
	done
done

mkdir /users/mturchin/data/1000G/mturchin20/Analyses/CAFE/Vs1/SubFiles/PermFiles

for i in `cat <(echo "Loh2017;3169 Neale2017;9727 GIANT2014_5;27673" | perl -lane 'print join("\n", @F);') | head -n 1 | tail -n 1`; do
	iVal1=`echo $i | perl -ane 'my @vals1 = split(/;/, $F[0]); print $vals1[0];'`; iSeed1=`echo $i | perl -ane 'my @vals1 = split(/;/, $F[0]); print $vals1[1];'`;
	for j in `cat <(echo "Height;5638" | perl -lane 'print join("\n", @F);') | head -n 1 | tail -n 1`; do
		jVal1=`echo $j | perl -ane 'my @vals1 = split(/;/, $F[0]); print $vals1[0];'`; jSeed1=`echo $j | perl -ane 'my @vals1 = split(/;/, $F[0]); print $vals1[1];'`;
		for k in `cat <(echo "lt5eNeg9;2759 lt1eNeg4;78364 lt5eNeg8;3869" | perl -lane 'print join("\n", @F);') | head -n 1 | tail -n 1`; do
			kVal1=`echo $k | perl -ane 'my @vals1 = split(/;/, $F[0]); print $vals1[0];'`; kSeed1=`echo $k | perl -ane 'my @vals1 = split(/;/, $F[0]); print $vals1[1];'`;
			TotalSeed1=$((iSeed1+jSeed1+kSeed1));
			echo $i $j $k $TotalSeed1

			if [ ! -d /users/mturchin/data/1000G/mturchin20/Analyses/CAFE/Vs1/SubFiles/PermFiles/$iVal1.$jVal1.$kVal1 ]; then
				mkdir /users/mturchin/data/1000G/mturchin20/Analyses/CAFE/Vs1/SubFiles/PermFiles/$iVal1.$jVal1.$kVal1
			fi

			rm -f /users/mturchin/data/1000G/mturchin20/Analyses/CAFE/Vs1/SubFiles/PermFiles/$iVal1.$jVal1.$kVal1/*
			zcat /users/mturchin/LabMisc/RamachandranLab/InterPath/CAFE/Data/GWASsnps/$iVal1.$jVal1.$kVal1.bmass.GrdyClmp.rsIDs.w1000GInfo.txt.gz | R -q -e "Data1 <- read.table(file('stdin'), header=F); \
			for (i in 1:10) { print(i); \
				AncAF <- round(as.numeric(as.character(Data1[i,8])), digits=2); Data2 <- c(); File3Count <- c(); for (j in c(-.01,0,.01)) { Filename3 <- paste(\"/users/mturchin/data/1000G/mturchin20/Analyses/CAFE/Vs1/SubFiles/PerAFFiles/CEUGBRTSIESNYRIESN.chrAll.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.SNPs.noEURfix.edit.wMeanInfo.permPrep.AncAF_\", as.character(AncAF + j), \".frq.gz\", sep=\"\"); if (file.exists(Filename3)) { Data3 <- read.table(Filename3, header=F); Data2 <- rbind(Data2, Data3); File3Count <- c(File3Count, j); }; }; print(c(AncAF, paste(File3Count, collapse=\",\"), dim(Data3))); \ 
				set.seed(as.numeric(as.character($TotalSeed1))+i); RowVals1 <- sample(1:nrow(Data3)); for (k in 1:10) { \
						write.table(Data3[RowVals1[k],], file=paste(\"/users/mturchin/data/1000G/mturchin20/Analyses/CAFE/Vs1/SubFiles/PermFiles/$iVal1.$jVal1.$kVal1/CEUGBRTSIESNYRIESN.chrAll.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.SNPs.noEURfix.edit.wMeanInfo.perm\", as.character(k), \".txt\", sep=\"\"), append=TRUE, quote=FALSE, row.names=FALSE, col.names=FALSE); \
				}; \
			};" > /users/mturchin/data/1000G/mturchin20/Analyses/CAFE/Vs1/SubFiles/PermFiles/$iVal1.$jVal1.$kVal1/CEUGBRTSIESNYRIESN.chrAll.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.SNPs.noEURfix.edit.wMeanInfo.permAll.Summary.txt
			gzip -f /users/mturchin/data/1000G/mturchin20/Analyses/CAFE/Vs1/SubFiles/PermFiles/$iVal1.$jVal1.$kVal1/*

		done
	done
done

for i in `cat <(echo "Height;1254 BMI;58923 Waist;49281 Hip;37485 WaistAdjBMI;82374 HipAdjBMI;6182" | perl -lane 'print join("\n", @F);') | grep -vE 'Waist;49|Hip;37' | tail -n 2 | head -n 1`; do
        for j in `cat <(echo $UKBioBankPops | perl -lane 'print join("\n", @F);') | grep -vE 'Ran10000|Irish' | head -n 4 | tail -n 1`; do
        for k in `cat <(echo "NonSyn Exonic ExonicPlus ExonicPlus20kb IntronicPlus20kb IntronicPlus20kb25 IntronicPlus20kb50 IntronicPlus20kb75 GD125000 GD500000 GD25000" | perl -lane 'print join("\n", @F);') | head -n 7 | tail -n 1`; do
                        ancestry1=`echo $j | perl -ane 'my @vals1 = split(/;/, $F[0]); print $vals1[0];'`; ancestry2=`echo $j | perl -ane 'my @vals1 = split(/;/, $F[0]); print $vals1[1];'`; NumSNPs=`zcat /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/Imputation/mturchin20/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.raw.edit.gz | head -n 1 | perl -ane 'print scalar(@F);'`; Pheno1=`echo $i | perl -ane 'my @vals1 = split(/;/, $F[0]); print $vals1[0];'`; PhenoSeed1=`echo $i | perl -ane 'my @vals1 = split(/;/, $F[0]); print $vals1[1];'`; AncSeed1=`echo $j | perl -ane 'my @vals1 = split(/;/, $F[0]); print $vals1[3];'`
#                       NumPaths=`cat /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.bim.AnnovarFormat.TableAnnovar.AAFix.hg19_multianno.GeneSNPs.SemiColonSplit.wRowPos.Regions.c2.${k}.noDups.txt | wc | awk '{ print $1 }'`
                        NumPaths=2592
                        echo $i $ancestry1 $ancestry2 $ancestry3 $k








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


~~~



