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

for pop1 in `echo "CEU GBR IBS TSI"`; do 
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

for chr1 in `echo {1..22} X`; do

	join <(zcat /users/mturchin/data/1000G/subsets/CEU/AFs/CEU.chr${chr1}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.SNPs.frq.gz | awk '{ print $2 "\t" $3 "\t" $4 "\t" $5}' | sort -k 1,1) <(zcat /users/mturchin/data/1000G/subsets/GBR/AFs/GBR.chr${chr1}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.SNPs.frq.gz | awk '{ print $2 "\t" $3 "\t" $4 "\t" $5}' | sort -k 1,1) | perl -lane 'if ($F[1] eq "0") { if ($F[4] eq "0") { if ($F[2] eq $F[5]) { print join("\t", @F); } } else { if ($F[2] eq $F[5]) { print join("\t", @F); } elsif ($F[2] eq $F[4]) { ($F[4], $F[5]) = ($F[5], $F[4]); $F[6] = 1 - $F[6]; print join("\t", @F); } else { $PH1 = 1; } } } else { if ($[4] eq "0") { if ($F[1] eq $F[5]) { ($F[4], $F[5]) = ($F[5], $F[4]); $F[6] = 1 - $F[6]; print join("\t", @F); } elsif ($F[2] eq $F[5]) { print join("\t", @F); } else { $PH1 = 1; } } else { if (($F[1] eq $F[4]) && ($F[2] eq $F[5])) { print join("\t", @F); } elsif (($F[1] eq $F[5]) && ($F[2] eq $F[4])) { ($F[4], $F[5]) = ($F[5], $F[4]); $F[6] = 1 - $F[6]; print join("\t", @F); } else { $PH1 = 1; } } }' 
	
	<(zcat /users/mturchin/data/1000G/subsets/CEU/AFs/CEU.chr${chr1}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.SNPs.frq.gz

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





