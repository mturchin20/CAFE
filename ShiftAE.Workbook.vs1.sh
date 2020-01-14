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

##zcat $HOME/data/mturchin/Neale2017/Vs2/50.assoc.tsv.gz | sed 's/:/ /g' | awk '{ print $5 "\t" $3 "\t" $4 "\t" $9 "\t" $10 "\t" $12 }' | grep -v pval | cat <(echo -e "rsid\ta1\ta2\tbeta\tse\tpval") - | gzip > $HOME/data/mturchin/Neale2017/Vs2/50.assoc.edits1.tsv.gz
zcat $HOME/data/mturchin/Neale2017/Vs2/50.assoc.tsv.gz | sed 's/:/ /g' | awk '{ print $5 "\t" $3 "\t" $4 "\t" $9 "\t" $10 "\t" $12 }' | grep -v pval | cat <(echo -e "rsid\ta1\ta2\tbeta\tse\tpval") - | gzip > $HOME/data/mturchin/Neale2017/Vs2/Neal2017.Height.edits.txt.gz
zcat $HOME/data/mturchin/Loh2017/body_HEIGHTz.sumstats.gz | awk '{ if ($4 != $6) { $8 = -1 * $8; val1 = $4; val2 = $5; $4 = val2; $5 = val1; }; print $1 "\t" $4 "\t" $5 "\t" $8 "\t" $9 "\t" $10 }' | grep -v SNP | cat <(echo -e "rsid\ta1\ta2\tbeta\tse\tpval") - | gzip > $HOME/data/mturchin/Loh2017/Loh2017.Height.edits.txt.gz
zcat /users/mturchin/Data2/GIANT/GIANT_HEIGHT_Wood_et_al_2014_publicrelease_HapMapCeuFreq.wUCSCGenomeBrowser_dbSNP130.vs1.txt.gz | awk '{ print $1 "\t" $2 "\t" $3 "\t" $5 "\t" $6 "\t" $7 }' | grep -v Allele | cat <(echo -e "rsid\ta1\ta2\tbeta\tse\tpval") - | gzip > /users/mturchin/Data2/GIANT/GIANT2014_5.Height.edits.txt.gz

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




PROBLEM FIX: 
[  mturchin@node1149  ~]$zcat /users/mturchin/data/1000G/mturchin20/Analyses/CAFE/Vs1/SubFiles/PerAFFiles/CEUGBRTSIESNYRIESN.chrAll.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.SNPs.noEURfix.edit.wMeanInfo.permPrep.AncAF_0.25.frq.gz | awk '{ print $1 }' | wc        
 117893  117893 1223855
[  mturchin@node1149  ~]$zcat /users/mturchin/data/1000G/mturchin20/Analyses/CAFE/Vs1/SubFiles/PerAFFiles/CEUGBRTSIESNYRIESN.chrAll.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.SNPs.noEURfix.edit.wMeanInfo.permPrep.AncAF_0.25.frq.gz | awk '{ print $1 }' | sort | uniq | wc
 117881  117881 1223720
[  mturchin@node1149  ~]$ls -lrt /users/mturchin/data/1000G/mturchin20/Analyses/CAFE/Vs1/SubFiles/PerAFFiles/CEUGBRTSIESNYRIESN.chrAll.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.SNPs.noEURfix.edit.wMeanInfo.permPrep.AncAF_0.69.frq.gz
-rw-r--r--. 1 mturchin sramacha 239 Jun 16  2019 /users/mturchin/data/1000G/mturchin20/Analyses/CAFE/Vs1/SubFiles/PerAFFiles/CEUGBRTSIESNYRIESN.chrAll.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.SNPs.noEURfix.edit.wMeanInfo.permPrep.AncAF_0.69.frq.gz
[  mturchin@node1149  ~]$zcat /users/mturchin/data/1000G/mturchin20/Analyses/CAFE/Vs1/SubFiles/PerAFFiles/CEUGBRTSIESNYRIESN.chrAll.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.SNPs.noEURfix.edit.wMeanInfo.permPrep.AncAF_0.69.frq.gz
rs12186596 G A 0.3485 1 1 0.4206 0.692275 0.343775 0.307725 0.307725 0.271675
rs12186596 G A 0.3485 1 1 0.4206 0.692275 0.343775 0.307725 0.307725 0.271675
rs12186596 G A 0.3485 1 1 0.4206 0.692275 0.343775 0.307725 0.307725 0.271675
rs12186596 G A 0.3485 1 1 0.4206 0.692275 0.343775 0.307725 0.307725 0.271675
rs12186596 G A 0.3485 1 0.4065 1 0.68875 0.34025 0.31125 0.28225 0.31125
rs12186596 G A 0.3485 1 0.4065 1 0.68875 0.34025 0.31125 0.28225 0.31125
rs12186596 G A 0.3485 1 0.4065 1 0.68875 0.34025 0.31125 0.28225 0.31125
rs12186596 G A 0.3485 1 0.4065 1 0.68875 0.34025 0.31125 0.28225 0.31125





mkdir /users/mturchin/data/1000G/mturchin20/Analyses/CAFE/Vs1/SubFiles/PerAFFiles/Loh2017
mkdir /users/mturchin/data/1000G/mturchin20/Analyses/CAFE/Vs1/SubFiles/PerAFFiles/Neale2017
mkdir /users/mturchin/data/1000G/mturchin20/Analyses/CAFE/Vs1/SubFiles/PerAFFiles/GIANT2014_5

CEUGBRTSIESNYRIESN.chrAll.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.SNPs.noEURfix.edit.wMeanInfo.permPrep.AncAF_0.48.frq.gz

#Loh2017
add in loh/neale/giant loop
for i in `ls -lrt /users/mturchin/data/1000G/mturchin20/Analyses/CAFE/Vs1/SubFiles/PerAFFiles/ | awk '{ print $9 }' | head -n 10`; do
	AFFile1="/users/mturchin/data/1000G/mturchin20/Analyses/CAFE/Vs1/SubFiles/PerAFFiles/${i}";
	AFFile1NewPre=`echo $AFFile1 | sed 's/.frq.gz//g'`; AFFile1New="${AFFile1NewPre}.Loh2017.frq";
	echo $AFFile1 $AFFile1New

	match height file allele to affile, switch beta if needed, then check if beta positive, and switch everything if not
	join <(zcat $AFFile1 | sort -k 1,1) <(zcat $HOME/data/mturchin/Loh2017/Loh2017.Height.edits.txt.gz | awk '{  

done



[  mturchin@node1149  ~]$zcat $HOME/data/mturchin/Neale2017/Vs2/50.assoc.tsv.gz | head -n 10
variant rsid    nCompleteSamples        AC      ytx     beta    se      tstat   pval
5:43888254:C:T  rs13184706      336474  1.23213e+04     3.83186e+02     -1.33660e-02    6.48300e-03     -2.06170e+00    3.92375e-02
5:43888493:C:T  rs58824264      336474  2.42483e+03     3.70381e+01     -1.88438e-02    1.44299e-02     -1.30589e+00    1.91592e-01
5:43888556:T:C  rs72762387      336474  1.64428e+04     6.65756e+02     -9.45691e-03    5.62698e-03     -1.68064e+00    9.28345e-02
5:43888648:C:T  rs115032754     336474  1.35047e+04     5.66843e+02     1.06178e-03     6.29484e-03     1.68674e-01     8.66053e-01
5:43888690:C:G  rs147555725     336474  1.24755e+03     4.51586e+01     -4.19957e-03    2.06522e-02     -2.03347e-01    8.38864e-01
5:43888838:G:C  rs13185925      336474  2.33424e+04     8.95665e+02     -8.69909e-03    4.74596e-03     -1.83294e+00    6.68117e-02
5:43889057:C:T  rs13189727      336474  1.25266e+04     4.21457e+02     -1.25947e-02    6.42476e-03     -1.96034e+00    4.99575e-02
5:43889207:A:G  rs4516856       336474  6.69113e+05     2.89791e+04     1.96508e-02     1.14359e-02     1.71835e+00     8.57338e-02
5:43889333:G:T  rs114787943     336474  3.03587e+03     5.56253e+01     -1.58425e-02    1.28856e-02     -1.22947e+00    2.18896e-01
[  mturchin@node1149  ~]$zcat /users/mturchin/Data2/GIANT/GIANT_HEIGHT_Wood_et_al_2014_publicrelease_HapMapCeuFreq.wUCSCGenomeBrowser_dbSNP130.vs1.txt.gz | head -n 10
MarkerName      Allele1 Allele2 Freq.Allele1.HapMapCEU  b       SE      p       N       ChrBP
rs4747841       A       G       0.551   -0.0011 0.0029  0.70    253213  10_10000134
rs4749917       T       C       0.436   0.0011  0.0029  0.70    253213  10_10000264
rs737656        A       G       0.367   -0.0062 0.0030  0.042   253116  10_100002728
rs737657        A       G       0.358   -0.0062 0.0030  0.041   252156  10_100002879
rs7086391       T       C       0.12    -0.0087 0.0038  0.024   248425  10_100003552
rs878177        T       C       0.3     0.014   0.0031  8.2e-06 251271  10_100003804
rs878178        A       T       0.644   0.0067  0.0031  0.029   253086  10_100003966
rs12219605      T       G       0.427   0.0011  0.0029  0.70    253213  10_10000458
rs3763688       C       G       0.144   -0.0022 0.0045  0.62    253056  10_100005552
[  mturchin@node1149  ~]$zcat $HOME/data/mturchin/Loh2017/body_HEIGHTz.sumstats.gz | head -n 10
SNP     CHR     POS     A1      A2      REF     EAF     Beta    se      P       N       INFO
rs10399793      1       49298   T       C       T       0.376205        -0.00104841     0.00272961      8.6E-01 458303  0.342797
rs2462492       1       54676   C       T       C       0.599409        -0.00505746     0.00270391      1.1E-01 458303  0.340158
rs3107975       1       55326   T       C       T       0.991557        0.00918983      0.015069        6.3E-01 458303  0.324228
rs74447903      1       57033   T       C       T       0.998221        -0.00814975     0.0335016       7.8E-01 458303  0.296256
1:70728_C_T     1       70728   C       T       C       0.997834        -0.0222717      0.0272397       5.9E-01 458303  0.365713
rs2462495       1       79033   A       G       A       0.00129168      -0.0429876      0.0346811       3.6E-01 458303  0.536566
rs114608975     1       86028   T       C       T       0.896401        -0.000955315    0.00432242      5.1E-01 458303  0.340885
rs6702460       1       91536   G       T       G       0.543017        -0.00371767     0.00266258      3.6E-01 458303  0.340746
rs8179466       1       234313  C       T       C       0.925496        0.00195976      0.00525025      7.0E-01 458303  0.311447



[  mturchin@node1149  ~]$zcat /users/mturchin/data/1000G/mturchin20/Analyses/CAFE/Vs1/CEUGBRTSIESNYRIESN.chr${chr1}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.SNPs.noEURfix.frq.gz | head -n 10
SNP     A1      A2      MAF     A1      A2      MAF     A1      A2      MAF     A1      A2      MAF     A1      A2      MAF     A1      A2      MAF
rs1000039       T       G       0.0202  T       G       0.03846 T       G       0.03738 T       G       0.02804 T       G       0.1898  T       G       0.1667
rs1000107       C       T       0.0303  C       T       0.01648 C       T       0.04206 C       T       0.04673 C       T       0.02315 C       T       0.01515
rs1000135       C       T       0.2273  C       T       0.1758  C       T       0.1636  C       T       0.1636  C       T       0.06019 C       T       0.07576
rs1000218       C       T       0.399   C       T       0.3901  C       T       0.3692  C       T       0.4299  C       T       0.3426  C       T       0.2879
rs1000219       T       G       0.399   T       G       0.3901  T       G       0.3692  T       G       0.4299  T       G       0.3565  T       G       0.298
rs1000251       T       C       0.1667  T       C       0.2363  T       C       0.229   T       C       0.1355  0       C       0       0       C       0
rs1000252       C       G       0.3333  C       G       0.4286  C       G       0.4206  C       G       0.3084  C       G       0.2685  C       G       0.298
rs1000266       C       G       0.2222  C       G       0.1813  C       G       0.1822  C       G       0.2196  C       G       0.3009  C       G       0.2828
rs1000278       T       C       0.4848  T       C       0.4396  T       C       0.3832  T       C       0.4112  T       C       0.2454  T       C       0.202
[  mturchin@node1149  ~]$zcat /users/mturchin/data/1000G/mturchin20/Analyses/CAFE/Vs1/CEUGBRTSIESNYRIESN.chr${chr1}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.SNPs.noEURfix.edit.wMeanInfo.frq.gz | head -n 10
rs1000039 T G 0.0202 0.03846 0.03738 0.02804 0.03102 0.01082 0.00744 0.00636 0.00298
rs1000107 C T 0.0303 0.01648 0.04206 0.04673 0.0338925 0.0035925 0.0174125 0.0081675 0.0128375
rs1000135 C T 0.2273 0.1758 0.1636 0.1636 0.182575 0.044725 0.00677499999999998 0.018975 0.018975
rs1000218 C T 0.399 0.3901 0.3692 0.4299 0.39705 0.00195000000000001 0.00695000000000001 0.02785 0.03285
rs1000219 T G 0.399 0.3901 0.3692 0.4299 0.39705 0.00195000000000001 0.00695000000000001 0.02785 0.03285
rs1000251 T C 0.1667 0.2363 0.229 0.1355 0.191875 0.025175 0.044425 0.037125 0.056375
rs1000252 C G 0.3333 0.4286 0.4206 0.3084 0.372725 0.039425 0.055875 0.047875 0.064325
rs1000266 C G 0.2222 0.1813 0.1822 0.2196 0.201325 0.020875 0.020025 0.019125 0.018275
rs1000278 T C 0.4848 0.4396 0.3832 0.4112 0.4297 0.0551 0.00990000000000002 0.0465 0.0185
rs1000280 A G 0.3485 0.3352 0.2757 0.229 0.2971 0.0514 0.0381 0.0214 0.0681
[  mturchin@node1149  ~]$zcat /users/mturchin/data/1000G/mturchin20/Analyses/CAFE/Vs1/SubFiles/PerAFFiles/CEUGBRTSIESNYRIESN.chrAll.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.SNPs.noEURfix.edit.wMeanInfo.permPrep.AncAF_0.53
.frq.gz | head -n 10
rs1006610 T C 0.4697 0.544 0.514 0.5888 0.529125 0.0594250000000001 0.014875 0.0151250000000001 0.0596749999999999
rs1007788 G C 0.4747 0.5165 0.5841 0.5374 0.528175 0.0534749999999999 0.011675 0.055925 0.00922500000000004
rs1007993 A C 0.4747 0.5165 0.5841 0.5374 0.528175 0.0534749999999999 0.011675 0.055925 0.00922500000000004
rs1014985 T G 0.4949 0.5495 0.5607 0.4953 0.5251 0.0302 0.0244 0.0356 0.0298
rs10157140 T C 0.4899 0.6044 0.5 0.5421 0.5341 0.0442 0.0703 0.0341 0.00800000000000001
rs10158342 A C 0.4798 0.511 0.5514 0.5654 0.5269 0.0471 0.0159 0.0245 0.0385
rs10158471 A G 0.4949 0.5055 0.5467 0.5654 0.528125 0.0332249999999999 0.022625 0.018575 0.0372750000000001
rs1032424 G A 0.4949 0.5275 0.5047 0.5981 0.5313 0.0364 0.00380000000000003 0.0266 0.0668
rs1038321 A G 0.4697 0.4945 0.6121 0.5561 0.5331 0.0634 0.0386 0.079 0.023
rs1041902 T A 0.4848 0.4835 0.5561 0.5841 0.527125 0.0423249999999999 0.043625 0.0289750000000001 0.056975
















mkdir /users/mturchin/LabMisc/RamachandranLab/InterPath/CAFE/Data
mkdir /users/mturchin/LabMisc/RamachandranLab/InterPath/CAFE/Data/GWASsnps

ln -s /users/mturchin/data/mturchin/Loh2017/body_HEIGHTz.sumstats.lt5eNeg9.bmass.GrdyClmp.rsIDs.gz /users/mturchin/LabMisc/RamachandranLab/InterPath/CAFE/Data/GWASsnps/Loh2017.Height.lt5eNeg9.bmass.GrdyClmp.rsIDs.gz
ln -s /users/mturchin/data/mturchin/Loh2017/body_HEIGHTz.sumstats.lt1eNeg4.bmass.GrdyClmp.rsIDs.gz /users/mturchin/LabMisc/RamachandranLab/InterPath/CAFE/Data/GWASsnps/Loh2017.Height.lt1eNeg4.bmass.GrdyClmp.rsIDs.gz
ln -s /users/mturchin/data/mturchin//Neale2017/Vs2/50.assoc.lt5eNeg9.bmass.GrdyClmp.rsIDs.gz /users/mturchin/LabMisc/RamachandranLab/InterPath/CAFE/Data/GWASsnps/Neale2017.Height.lt5eNeg9.bmass.GrdyClmp.rsIDs.gz
ln -s /users/mturchin/data/mturchin//Neale2017/Vs2/50.assoc.lt1eNeg4.bmass.GrdyClmp.rsIDs.gz /users/mturchin/LabMisc/RamachandranLab/InterPath/CAFE/Data/GWASsnps/Neale2017.Height.lt1eNeg4.bmass.GrdyClmp.rsIDs.gz
ln -s /users/mturchin/Data2/GIANT/GIANT_HEIGHT_Wood_et_al_2014_publicrelease_HapMapCeuFreq.wUCSCGenomeBrowser_dbSNP130.vs1.lt5eNeg9.GrdyClmp.rsIDs.gz /users/mturchin/LabMisc/RamachandranLab/InterPath/CAFE/Data/GWASsnps/GIANT2014_5.Height.lt5eNeg9.bmass.GrdyClmp.rsIDs.gz
ln -s /users/mturchin/Data2/GIANT/GIANT_HEIGHT_Wood_et_al_2014_publicrelease_HapMapCeuFreq.wUCSCGenomeBrowser_dbSNP130.vs1.lt5eNeg8.GrdyClmp.rsIDs.gz /users/mturchin/LabMisc/RamachandranLab/InterPath/CAFE/Data/GWASsnps/GIANT2014_5.Height.lt5eNeg8.bmass.GrdyClmp.rsIDs.gz
ln -s /users/mturchin/Data2/GIANT/GIANT_HEIGHT_Wood_et_al_2014_publicrelease_HapMapCeuFreq.wUCSCGenomeBrowser_dbSNP130.vs1.lt1eNeg4.GrdyClmp.rsIDs.gz /users/mturchin/LabMisc/RamachandranLab/InterPath/CAFE/Data/GWASsnps/GIANT2014_5.Height.lt1eNeg4.bmass.GrdyClmp.rsIDs.gz
	
for i in `cat <(echo "Loh2017 Neale2017 GIANT2014_5" | perl -lane 'print join("\n", @F);') | head -n 3 | tail -n 3`; do
	for j in `cat <(echo "Height" | perl -lane 'print join("\n", @F);') | head -n 1 | tail -n 1`; do
		for k in `cat <(echo "lt5eNeg9 lt1eNeg4 lt5eNeg8" | perl -lane 'print join("\n", @F);') | head -n 2 | tail -n 2`; do
			echo $i $j $k

#			join <(zcat /users/mturchin/LabMisc/RamachandranLab/InterPath/CAFE/Data/GWASsnps/$i.$j.$k.bmass.GrdyClmp.rsIDs.gz | sort) <(zcat /users/mturchin/data/1000G/mturchin20/Analyses/CAFE/Vs1/CEUGBRTSIESNYRIESN.chr*.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.SNPs.noEURfix.edit.wMeanInfo.frq.gz | sort -k 1,1) | gzip > /users/mturchin/LabMisc/RamachandranLab/InterPath/CAFE/Data/GWASsnps/$i.$j.$k.bmass.GrdyClmp.rsIDs.w1000GInfo.txt.gz
			zcat /users/mturchin/LabMisc/RamachandranLab/InterPath/CAFE/Data/GWASsnps/$i.$j.$k.bmass.GrdyClmp.rsIDs.w1000GInfo.txt.gz | R -q -e "Data1 <- read.table(file('stdin'), header=F); write.table(matrix(c(sum(abs(mean(Data1[,ncol(Data1)-3]) - Data1[,ncol(Data1)-3])), sum(abs(mean(Data1[,ncol(Data1)-2]) - Data1[,ncol(Data1)-2])), sum(abs(mean(Data1[,ncol(Data1)-1]) - Data1[,ncol(Data1)-1])), sum(abs(mean(Data1[,ncol(Data1)]) - Data1[,ncol(Data1)]))), nrow=1), quote=FALSE, col.names=FALSE, row.names=FALSE);"	

		done
	done
done

mkdir /users/mturchin/data/1000G/mturchin20/Analyses/CAFE/Vs1/SubFiles/PermFiles

#From: https://stackoverflow.com/questions/14904983/how-do-i-check-the-existence-of-a-local-file
for i in `cat <(echo "Loh2017;3169 Neale2017;9727 GIANT2014_5;27673" | perl -lane 'print join("\n", @F);') | head -n 3 | tail -n 1`; do
	iVal1=`echo $i | perl -ane 'my @vals1 = split(/;/, $F[0]); print $vals1[0];'`; iSeed1=`echo $i | perl -ane 'my @vals1 = split(/;/, $F[0]); print $vals1[1];'`;
	for j in `cat <(echo "Height;5638" | perl -lane 'print join("\n", @F);') | head -n 1 | tail -n 1`; do
		jVal1=`echo $j | perl -ane 'my @vals1 = split(/;/, $F[0]); print $vals1[0];'`; jSeed1=`echo $j | perl -ane 'my @vals1 = split(/;/, $F[0]); print $vals1[1];'`;
		for k in `cat <(echo "lt5eNeg9;2759 lt1eNeg4;78364 lt5eNeg8;3869" | perl -lane 'print join("\n", @F);') | head -n 2 | tail -n 2`; do
			kVal1=`echo $k | perl -ane 'my @vals1 = split(/;/, $F[0]); print $vals1[0];'`; kSeed1=`echo $k | perl -ane 'my @vals1 = split(/;/, $F[0]); print $vals1[1];'`;
			TotalSeed1=$((iSeed1+jSeed1+kSeed1));
			echo $i $j $k $TotalSeed1

			if [ ! -d /users/mturchin/data/1000G/mturchin20/Analyses/CAFE/Vs1/SubFiles/PermFiles/$iVal1.$jVal1.$kVal1 ]; then
				mkdir /users/mturchin/data/1000G/mturchin20/Analyses/CAFE/Vs1/SubFiles/PermFiles/$iVal1.$jVal1.$kVal1
			fi

			rm -f /users/mturchin/data/1000G/mturchin20/Analyses/CAFE/Vs1/SubFiles/PermFiles/$iVal1.$jVal1.$kVal1/*
			zcat /users/mturchin/LabMisc/RamachandranLab/InterPath/CAFE/Data/GWASsnps/$iVal1.$jVal1.$kVal1.bmass.GrdyClmp.rsIDs.w1000GInfo.txt.gz | R -q -e "Data1 <- read.table(file('stdin'), header=F); \
			for (i in 1:nrow(Data1)) { \
				AncAF <- round(as.numeric(as.character(Data1[i,8])), digits=2); Data2 <- c(); File3Count <- c(); for (j in c(-.01,0,.01)) { \
					Filename3 <- paste(\"/users/mturchin/data/1000G/mturchin20/Analyses/CAFE/Vs1/SubFiles/PerAFFiles/CEUGBRTSIESNYRIESN.chrAll.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.SNPs.noEURfix.edit.wMeanInfo.permPrep.AncAF_\", as.character(AncAF + j), \".frq.gz\", sep=\"\"); if (file.exists(Filename3)) { Data3 <- read.table(Filename3, header=F); Data2 <- rbind(Data2, Data3); File3Count <- c(File3Count, j); }; }; \
				cat(c(i, as.numeric(as.character($TotalSeed1))+i, AncAF, paste(File3Count, collapse=\",\"), dim(Data2)), \"\n\", sep=\"\t\"); \ 
				set.seed(as.numeric(as.character($TotalSeed1))+i); RowVals1 <- sample(1:nrow(Data2)); for (k in 1:1000) { \
						write.table(Data2[RowVals1[k],], file=paste(\"/users/mturchin/data/1000G/mturchin20/Analyses/CAFE/Vs1/SubFiles/PermFiles/$iVal1.$jVal1.$kVal1/CEUGBRTSIESNYRIESN.chrAll.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.SNPs.noEURfix.edit.wMeanInfo.perm\", as.character(k), \".txt\", sep=\"\"), append=TRUE, quote=FALSE, row.names=FALSE, col.names=FALSE); \
				}; \
			};" | grep -v ^\> > /users/mturchin/data/1000G/mturchin20/Analyses/CAFE/Vs1/SubFiles/PermFiles/$iVal1.$jVal1.$kVal1/CEUGBRTSIESNYRIESN.chrAll.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.SNPs.noEURfix.edit.wMeanInfo.permAll.Summary.txt
			gzip -f /users/mturchin/data/1000G/mturchin20/Analyses/CAFE/Vs1/SubFiles/PermFiles/$iVal1.$jVal1.$kVal1/*

		done
	done
done

##mkdir /users/mturchin/data/1000G/mturchin20/Analyses/CAFE/Vs1/Analyses
##mkdir /users/mturchin/data/1000G/mturchin20/Analyses/CAFE/Vs1/Analyses/Perms
mkdir /users/mturchin/LabMisc/RamachandranLab/InterPath/CAFE/Vs1/Analyses/Perms

for i in `cat <(echo "Loh2017;3169 Neale2017;9727 GIANT2014_5;27673" | perl -lane 'print join("\n", @F);') | head -n 3 | tail -n 3`; do
        iVal1=`echo $i | perl -ane 'my @vals1 = split(/;/, $F[0]); print $vals1[0];'`; iSeed1=`echo $i | perl -ane 'my @vals1 = split(/;/, $F[0]); print $vals1[1];'`;
        for j in `cat <(echo "Height;5638" | perl -lane 'print join("\n", @F);') | head -n 1 | tail -n 1`; do
                jVal1=`echo $j | perl -ane 'my @vals1 = split(/;/, $F[0]); print $vals1[0];'`; jSeed1=`echo $j | perl -ane 'my @vals1 = split(/;/, $F[0]); print $vals1[1];'`;
                for k in `cat <(echo "lt5eNeg9;2759 lt1eNeg4;78364 lt5eNeg8;3869" | perl -lane 'print join("\n", @F);') | head -n 2 | tail -n 2`; do
                        kVal1=`echo $k | perl -ane 'my @vals1 = split(/;/, $F[0]); print $vals1[0];'`; kSeed1=`echo $k | perl -ane 'my @vals1 = split(/;/, $F[0]); print $vals1[1];'`;
                        TotalSeed1=$((iSeed1+jSeed1+kSeed1));
                        echo $i $j $k $TotalSeed1

			rm -f /users/mturchin/LabMisc/RamachandranLab/InterPath/CAFE/Vs1/Analyses/Perms/$iVal1.$jVal1.$kVal1.CEUGBRTSIESNYRIESN.chrAll.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.SNPs.noEURfix.edit.wMeanInfo.permAll.Results.txt; for l in {1..20}; do
				zcat /users/mturchin/data/1000G/mturchin20/Analyses/CAFE/Vs1/SubFiles/PermFiles/$iVal1.$jVal1.$kVal1/CEUGBRTSIESNYRIESN.chrAll.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.SNPs.noEURfix.edit.wMeanInfo.perm${l}.txt.gz | R -q -e "Data1 <- read.table(file('stdin'), header=F); write.table(matrix(c(sum(abs(mean(Data1[,ncol(Data1)-3]) - Data1[,ncol(Data1)-3])), sum(abs(mean(Data1[,ncol(Data1)-2]) - Data1[,ncol(Data1)-2])), sum(abs(mean(Data1[,ncol(Data1)-1]) - Data1[,ncol(Data1)-1])), sum(abs(mean(Data1[,ncol(Data1)]) - Data1[,ncol(Data1)]))), nrow=1), quote=FALSE, col.names=FALSE, row.names=FALSE);" | grep -v ^\> >> /users/mturchin/LabMisc/RamachandranLab/InterPath/CAFE/Vs1/Analyses/Perms/$iVal1.$jVal1.$kVal1.CEUGBRTSIESNYRIESN.chrAll.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.SNPs.noEURfix.edit.wMeanInfo.permAll.Results.txt
			done 

                done
        done
done

#				...Data1[,ncol(Data1)-3]), sd(Data1[,ncol(Data1)-2]), sd(Data1[,ncol(Data1)-1]), sd(Data1[,ncol(Data1)])), nrow=1), quote=FALSE, col.names=FALSE, row.names=FALSE);" | grep -v ^\> >> /users/mturchin/LabMisc/RamachandranLab/InterPath/CAFE/Vs1/Analyses/Perms/$iVal1.$jVal1.$kVal1.CEUGBRTSIESNYRIESN.chrAll.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.SNPs.noEURfix.edit.wMeanInfo.permAll.Results.txt 

for i in `cat <(echo "Loh2017;3169 Neale2017;9727 GIANT2014_5;27673" | perl -lane 'print join("\n", @F);') | head -n 3 | tail -n 3`; do
        iVal1=`echo $i | perl -ane 'my @vals1 = split(/;/, $F[0]); print $vals1[0];'`; iSeed1=`echo $i | perl -ane 'my @vals1 = split(/;/, $F[0]); print $vals1[1];'`;
        for j in `cat <(echo "Height;5638" | perl -lane 'print join("\n", @F);') | head -n 1 | tail -n 1`; do
                jVal1=`echo $j | perl -ane 'my @vals1 = split(/;/, $F[0]); print $vals1[0];'`; jSeed1=`echo $j | perl -ane 'my @vals1 = split(/;/, $F[0]); print $vals1[1];'`;
                for k in `cat <(echo "lt5eNeg9;2759 lt1eNeg4;78364 lt5eNeg8;3869" | perl -lane 'print join("\n", @F);') | head -n 2 | tail -n 2`; do
                        kVal1=`echo $k | perl -ane 'my @vals1 = split(/;/, $F[0]); print $vals1[0];'`; kSeed1=`echo $k | perl -ane 'my @vals1 = split(/;/, $F[0]); print $vals1[1];'`;
                        TotalSeed1=$((iSeed1+jSeed1+kSeed1));
                        echo $i $j $k $TotalSeed1

			cat /users/mturchin/LabMisc/RamachandranLab/InterPath/CAFE/Vs1/Analyses/Perms/$iVal1.$jVal1.$kVal1.CEUGBRTSIESNYRIESN.chrAll.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.SNPs.noEURfix.edit.wMeanInfo.permAll.Results.txt | head -n 10

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



~~~


