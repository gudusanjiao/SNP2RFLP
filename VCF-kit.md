# Workflow of using VCF-kit to design primers related to polymorphic restriction site in vcf files
Nan Hu / 
Nov. 22, 2020

---

## Installation and operating environments for VCF-kit (important!)
To prepare for installing VCF-kit, make sure system has pip and anaconda installed ahead. 
If new users do not have these two tools installed, you can refer to: [pip installing guide](https://pip.pypa.io/en/stable/installing/) 
and [anaconda installing guide](https://docs.anaconda.com/anaconda/install/linux/) to get instructions for installation. Once required tools are prepared, 
execute below commands in commandline:
''
conda config --add channels bioconda
conda config --add channels conda-forge
conda create -n vcf-kit \
  danielecook::vcf-kit=0.2.6 \
  "bwa>=0.7.17" \
  "samtools>=1.10" \
  "bcftools>=1.10" \
  "blast>=2.2.31" \
  "muscle>=3.8.31" \
  "primer3>=2.5.0"
''
