# Workflow of using VCF-kit to design primers related to polymorphic restriction site in vcf files
Nan Hu / 
Nov. 22, 2020

---

## Installation and operating environments for VCF-kit (important!)
To prepare for installing VCF-kit, make sure system has pip and anaconda installed ahead. 
If new users do not have these two tools installed, you can refer to: [pip installing guide](https://pip.pypa.io/en/stable/installing/) 
and [anaconda installing guide](https://docs.anaconda.com/anaconda/install/linux/) to get instructions for installation. Once required tools are prepared, 
execute below commands in commandline:
```bash
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
```
These process a virtual enviroment only for VCF-kit. We can activate this enviroment by:
```bash
conda activate vcf-kit
```
Since this software have not been updated for several years, there are some python codes not accommodating new biopython package. 
We need to install old versions of biopython and fix some issues there.
```bash
# uninstall new version of biopython
conda uninstall biopython
# install v1.77 biopython
conda install biopython=1.77
# uninstall numpy thoroughtly, repeat this step until there is nothing to uninstall
conda uninstall numpy
# re-install numpy
conda install numpy
```
Check the availability of VCF-kit with our expected function by:
```bash
vk primer --help
```
Last step to prepare the environment is to fix some coding issue in finding primer3 software. 
VCF-kit expected we installed Primer3 under system folder so by default it will search primer3-config files in some of these folders.
We need to find where our primer3-config files are and manually add to the searching range.
```bash
find ~ -iname primer3_config -type d
```
There will be a result at ```~/anaconda3/lib/python3.6/site-packages/vcfkit/static/primer3_config``` Save this path and do another search:
```bash
find ~ -iname primer3.py -type f
```
Use text editors (like: nano, vi etc.) to open the file locate at ```/home/nan/anaconda3/envs/vcf-kit/lib/python3.7/site-packages/vcfkit/utils/primer3.py```
In this python script, find 'class primer3'. Under \# Global default, there is a block of script:
```python
        # Global default
        thermo_paths = ["/usr/local/share/primer3_config/",
                        "/usr/local/share/primer3/primer3_config/",
                        "~/.linuxbrew/share/primer3_config/",
                        "~/.linuxbrew/share/primer3/primer3_config/",
                        "/.linuxbrew/share/primer3_config/",
                        "/.linuxbrew/share/primer3/primer3_config/",
                        primer3_config]
```
Add our saved path to this constant default. Completed one will look like:
```python
        # Global default
        thermo_paths = ["/usr/local/share/primer3_config/",
                        "/usr/local/share/primer3/primer3_config/",
                        "~/.linuxbrew/share/primer3_config/",
                        "~/.linuxbrew/share/primer3/primer3_config/",
                        "/.linuxbrew/share/primer3_config/",
                        "/.linuxbrew/share/primer3/primer3_config/",
                        "/home/nan/anaconda3/lib/python3.6/site-packages/vcfkit/static/primer3_config",
                        primer3_config]
```
Now, we have finished enviroment settings for VCF-kit primer function.

## 






