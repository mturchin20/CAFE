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

mkdir /users/mturchin/LabMisc/RamachandranLab/InterPath/CAFE
mkdir /users/mturchin/LabMisc/RamachandranLab/InterPath/CAFE/Vs1
mkdir /users/mturchin/LabMisc/RamachandranLab/InterPath/CAFE/Vs1/Analyses

for pop1 in `echo "CEU GBR TSI IBS"`; do 
	echo $pop1

	if [ ! -d /users/mturchin/data/1000G/subsets/$pop1/AFs ] ; then
		mkdir /users/mturchin/data/1000G/subsets/$pop1/AFs
	fi
	
	for chr1 in `echo {1..22} X`; do

#		plink --bfile /users/mturchin/data/1000G/subsets/$pop1/${pop1}.chr${chr1}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.SNPs --freq --out /users/mturchin/data/1000G/subsets/$pop1/AFs/${pop1}.chr${chr1}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.SNPs  
		gzip -f /users/mturchin/data/1000G/subsets/$pop1/AFs/${pop1}.chr${chr1}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.SNPs.frq

	done
done

mkdir /users/mturchin/data/1000G/mturchin20/Analyses
mkdir /users/mturchin/data/1000G/mturchin20/Analyses/CAFE
mkdir /users/mturchin/data/1000G/mturchin20/Analyses/CAFE/Vs1

#NOTE -- all 1000G files are on the same strand, and since pop subsets were extracted from the main source AllPops files, they should all start with the same set of SNPs at least (read: http://www.internationalgenome.org/faq/what-strand-are-variants-your-vcf-file/)

for chr1 in `echo {1..1}`; do

	join <(zcat /users/mturchin/data/1000G/subsets/CEU/AFs/CEU.chr${chr1}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.SNPs.frq.gz | awk '{ if ($6 > 125) { print $0 } } ' | awk '{ print $2 "\t" $3 "\t" $4 "\t" $5}' | sort -k 1,1) <(zcat /users/mturchin/data/1000G/subsets/GBR/AFs/GBR.chr${chr1}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.SNPs.frq.gz | awk '{ if ($6 > 125) { print $0 } } ' | awk '{ print $2 "\t" $3 "\t" $4 "\t" $5}' | sort -k 1,1) | perl -lane 'if ($F[1] eq "0") { if ($F[4] eq "0") { if ($F[2] eq $F[5]) { print join("\t", @F); } } else { if ($F[2] eq $F[5]) { print join("\t", @F); } elsif ($F[2] eq $F[4]) { ($F[4], $F[5]) = ($F[5], $F[4]); $F[6] = 1 - $F[6]; print join("\t", @F); } else { $PH1 = 1; } } } else { if ($F[4] eq "0") { if ($F[1] eq $F[5]) { ($F[4], $F[5]) = ($F[5], $F[4]); $F[6] = 1 - $F[6]; print join("\t", @F); } elsif ($F[2] eq $F[5]) { print join("\t", @F); } else { $PH1 = 1; } } else { if (($F[1] eq $F[4]) && ($F[2] eq $F[5])) { print join("\t", @F); } elsif (($F[1] eq $F[5]) && ($F[2] eq $F[4])) { ($F[4], $F[5]) = ($F[5], $F[4]); $F[6] = 1 - $F[6]; print join("\t", @F); } else { $PH1 = 1; } } }' | \ 
	join - <(zcat /users/mturchin/data/1000G/subsets/TSI/AFs/TSI.chr${chr1}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.SNPs.frq.gz | awk '{ if ($6 > 125) { print $0 } } ' | awk '{ print $2 "\t" $3 "\t" $4 "\t" $5}' | sort -k 1,1) | perl -lane 'if ($F[1] eq "0") { if ($F[7] eq "0") { if ($F[2] eq $F[8]) { print join("\t", @F); } } else { if ($F[2] eq $F[8]) { print join("\t", @F); } elsif ($F[2] eq $F[7]) { ($F[7], $F[8]) = ($F[8], $F[7]); $F[9] = 1 - $F[9]; print join("\t", @F); } else { $PH1 = 1; } } } else { if ($F[7] eq "0") { if ($F[1] eq $F[8]) { ($F[7], $F[8]) = ($F[8], $F[7]); $F[9] = 1 - $F[9]; print join("\t", @F); } elsif ($F[2] eq $F[8]) { print join("\t", @F); } else { $PH1 = 1; } } else { if (($F[1] eq $F[7]) && ($F[2] eq $F[8])) { print join("\t", @F); } elsif (($F[1] eq $F[8]) && ($F[2] eq $F[7])) { ($F[7], $F[8]) = ($F[8], $F[7]); $F[9] = 1 - $F[9]; print join("\t", @F); } else { $PH1 = 1; } } }' | \
	
	head -n 10 
	\ 

	4 -> 7
	5 -> 8
	6 -> 9

	SNP A1 A2 MAF A1 A2 MAF A1 A2 MAF
	rs1000033 G T 0.2071 G T 0.2033 G T 0.1682
	rs1000050 C T 0.1364 C T 0.1429 C T 0.1449
	rs1000070 T C 0.2828 T C 0.3626 T C 0.3224
	rs1000072 0 G 0 0 G 0 0 G 0

	
	join - <(zcat /users/mturchin/data/1000G/subsets/IBS/AFs/IBS.chr${chr1}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.SNPs.frq.gz | awk '{ print $2 "\t" $3 "\t" $4 "\t" $5}' | sort -k 1,1) | head -n 10

done















mkdir /users/mturchin/LabMisc/RamachandranLab/InterPath
mkdir /users/mturchin/LabMisc/RamachandranLab/InterPath/Vs1
mkdir /users/mturchin/LabMisc/RamachandranLab/InterPath/Vs1/Analyses
cd /users/mturchin/LabMisc/RamachandranLab/InterPath; git clone https://github.com/mturchin20/Ramachandran_InterPath
mkdir /users/mturchin/data/mturchin/Broad
mkdir /users/mturchin/data/mturchin/Broad/MSigDB
#NOTE -- need to log-in with an e-mail address and then DL directly from webpage, so cannot really do wget setup (I think anyways); so files were downloaded to MacBook Air and then scp'd over
#cd /users/mturchin/data/mturchin/Broad/MSigDB
#wget
#From MackBook Air
#scp -p /Users/mturchin20/Downloads/*gmt mturchin@ssh.ccv.brown.edu:/users/mturchin/data/mturchin/Broad/MSigDB/.

#NOTE -- copy and pasted 'InterPath_Simulation.R' and 'InterPath.cpp' from Slack via Lorin

cp -p /users/mturchin/LabMisc/RamachandranLab/InterPath/InterPath.cpp /users/mturchin/LabMisc/RamachandranLab/InterPath/Vs1/InterPath.mtEdits.vs1.cpp
cp -p /users/mturchin/LabMisc/RamachandranLab/InterPath/InterPath_Simulation.R /users/mturchin/LabMisc/RamachandranLab/InterPath/Vs1/InterPath.mtEdits.vs1.R
cp -p /users/mturchin/LabMisc/RamachandranLab/InterPath/InterPath_Simulation.R /users/mturchin/LabMisc/RamachandranLab/InterPath/Vs1/InterPath.Source.vs1.R

#See /users/mturchin/.bash_profile for how eventually got the below g++ call to work/function to get the stand-alone armadillo c++ code going
#g++ Test1.cpp -o Test1.out  -O2 /gpfs/data/sramacha/mturchin/miniconda2/envs/FORCE/lib/libgfortran.so.4 /gpfs/data/sramacha/mturchin/miniconda2/envs/FORCE/lib/libarmadillo.so.8 -larmadillo -lopenblas -larpack -lgfortran -L/gpfs/data/sramacha/mturchin/miniconda2/envs/FORCE/lib
g++ Test1.cpp -o Test1.out  -O2 /gpfs/data/sramacha/mturchin/miniconda2/envs/InterPath/lib/libgfortran.so.4 /gpfs/data/sramacha/mturchin/miniconda2/envs/InterPath/lib/libarmadillo.so.8 -larmadillo -lopenblas -larpack -lgfortran -L/gpfs/data/sramacha/mturchin/miniconda2/envs/InterPath/lib

mkdir /users/mturchin/data/mturchin/20180611_annovar
mkdir /users/mturchin/data/mturchin/20180611_annovar/humandb





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


~~~



