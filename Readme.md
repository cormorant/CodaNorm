CodaNorm
**a software package for the body-wave attenuation calculation by coda-normalization method**
===

### Abstract

The software package CodaNorm is an open source seismological software. It consists of two separate utilities: a converter baikal2gse allows the conversion of the seismic waveforms of the regional seismic station networks in Southern Siberia and Yakutia to the GSE2 format and the program CodaNorm, designed to estimate the seismic wave attenuation of local and regional earthquakes and explosions using the extended coda-normalization method for different frequency ranges.
CodaNorm program allows to calculate the seismic quality factor using the coda-normalization method [Aki, 1980]. 

### Motivation and significance

The software package CodaNorm is an open source seismological software and it is designed to estimate the body-wave attenuation of local and regional earthquakes and explosions. An interpreted programming language Python was used to develop the software package.

The program CodaNorm allows the estimation of the seismic quality factor (Qp, Qs), its frequency dependence (n) and attenuation decrement () for body P- and S-waves by the coda-normalization method [Aki, 1980] for different frequency ranges selected by a user. The initial data used are: 1. The earthquake (or explosion) waveforms (seismograms) in the GSE2 format; 2. the text file containing the information on the earthquake origin time, coordinates, arrival times of direct body-waves (P- and S-waves), and the path to the waveform database as well. The information on the regional attenuation parameters (the seismic quality factor, its frequency dependence and attenuation decrement) is necessary to correct the decay shake model from the earthquakes on the traces from the seismically active zones in the main urban areas, as well as for the further calculation of synthetic accelerograms and the evaluation of the parameters of the vibration for the possible strong earthquakes.

### System requirements

Python installed to run the program (version 2.6 or more) is required [http://python.org].
Python 3 with Obspy installed (since version 0.10.1) is supported.
Installing Obspy implies Numpy, Scipy and Matplotlib libraries. For details on installing Obspy, please, consult the [official guide] (https://github.com/obspy/obspy/wiki]).

CodaNorm program compiled executables are available for Windows 32-bit. To start them Python with libraries installed is not required.
