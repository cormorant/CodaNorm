css2baikal
===

**A program to convert the data format CSS-3.0 into Байкал-5 (ХХ)**

### Usage

View all available options with the command

**css2baikal.exe -h**

or

**css2baikal.exe --help**

Run data conversion from CSS-3.0 in a specified folder:

**css2baikal.exe FOLDER_CSS**

for example:

**css2baikal.exe ORYE**

Will find all file for converting in folder ORYE and save the result in folder XX (by default).
To specify the output folder option is used **-o (--out) option**.

Other options:

**-m (--merge)**

Files for each station and the channel will be combined, the output files will break for N minutes (default 10 minutes, an option ** - n ** changes the length of the file).

Command:

**css2baikal.exe ORYE -m**

Will write files from all channels combined in the output folder + ORYE.

Name output files will be written in format **'Y-M-D-H-M-S.f'**.
