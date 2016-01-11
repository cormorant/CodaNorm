baikal2gse
===

**A program to convert the data format Baikal-5 (ХХ) into GSE2**

### Abstract

Waveforms of earthquakes recorded are stored in a regional format Baykal-5. This script converts data from Baykal-5 format to GSE2 format.



### Usage

View all of the options available with the command

**baikal2gse.exe -h**



baikal2gse.py [-h] [-V] [-o OUTDIR] dirs [dirs ...]

required arguments:
  dirs                  path to directory with waveform data (Baikal-5)

optional arguments:
  -h, --help            help
  -V, --version         version
  -o OUTDIR, --outdir OUTDIR
                        the path to save the output data (default is "gse2")
