CodaNorm
===

**a software package for the body-wave attenuation calculation by coda-normalization method**


### Abstract

The software package CodaNorm is an open source seismological software. It consists of two separate utilities: a converter baikal2gse allows the conversion of the seismic waveforms of the regional seismic station networks in Southern Siberia and Yakutia to the GSE2 format and the program CodaNorm, designed to estimate the seismic wave attenuation of local and regional earthquakes and explosions using the extended coda-normalization method for different frequency ranges.
CodaNorm program allows to calculate the seismic quality factor using the coda-normalization method [Aki, 1980]. 

### Motivation and significance

The software package CodaNorm is an open source seismological software and it is designed to estimate the body-wave attenuation of local and regional earthquakes and explosions. An interpreted programming language Python was used to develop the software package.

### System requirements

Python installed to run the program (version 2.6 or more) is required [http://python.org].
Python 3 with Obspy installed (since version 0.10.1) is supported.
Installing Obspy implies Numpy, Scipy and Matplotlib libraries. For details on installing Obspy, please, consult the [official guide] (https://github.com/obspy/obspy/wiki).

CodaNorm program compiled executables are available for Windows 32-bit. To start them Python with libraries installed is not required.

### References

* Aki, K. (1980). Attenuation of shear waves in the lithosphere for frequencies from 0.05 to 25 Hz, Phys. Earth Planet. Interiors, 21. Р. 50–60.
* Dobrynina A.A. Coda-wave attenuation in the Baikal rift system lithosphere // Phys. Earth Planet. In., 2011, v. 188, pp. 121-126. doi:10.1016/j.pepi.2011.05.008.
* M. Beyreuther, R. Barsch, L. Krischer, T. Megies, Y. Behr and J. Wassermann (2010). ObsPy: A Python Toolbox for Seismology. SRL, 81(3), 530-533. DOI: 10.1785/gssrl.81.3.530
* Predein P.A., Tubanov Ts.A., German E. Study of seismic wave attenuation in the crust of the Baikal rift by coda-normalization method // Lithosphere structure and geodynamics: Articles XXVI-Russian Youth Conference (Irkutsk, April 20–25, 2015). Irkutsk: Earth’s crust Institute, 2015. P. 140–142 (in Russian).
* Rautian T.G., Khalturin V.I., 1978. The use of coda for determination of the earthquake source spectrum. Bulletin of the Seismological Society of America 68 (4), 923–948.
