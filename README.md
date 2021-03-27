# SPECFIT

Specfit is as a Python/C++/ROOT utility for fitting arbitrary cosmic ray energy spectrum (flux) measurements into a model functions using binned Poisson log likelihood method and displaying the results.  Constant and nonlinear energy scale corrections are supported.

* Author: Dmitri Ivanov <dmiivanov@gmail.com>
* Release: 09/2020


## Requirements:
- Linux OS
- GCC >= 4.6.3 
- CERN ROOT >= 5.34/36
- Python >= 2.7.3
- PyROOT properly installed ("import ROOT" line in Python must work)

### Preferred:
- Linux OS
- gcc >= 8.3.1
- c++11 standard
- Python >= 3.6.8
- CERN ROOT >= 6.20/04

## Install / Compile with make:
```bash
cd specfit
make cleanall
make
```

### Optionally, add environment 

```bash
source /full/path/to/specfit/specfit_env.sh 
```	
to your ${HOME}/.bashrc file if you're using BASH. 

### Install and Compile with CMake:

Alternatively, you can use cmake >= 3.9 to build specfit with ROOT >= 5.34.01, provided that ROOT was build with cmake.

```bash
cd specfit
mkdir obj
cd obj
cmake ..
make -j3
```
## Usage:

### Python3 (Preferred): 

```bash
./specfit.py --help | more 
```
shows a detailed manual.

### Python2 (Compatible with): 

```bash
python2 ./specfit.py --help | more
```


## Tested:

* Ubuntu 20.04.2 LTS, g++ (Ubuntu 9.3.0-17ubuntu1~20.04) 9.3.0, Python 3.8.5, ROOT 6.23/9
* Ubuntu 20.04.2 LTS, g++ (Ubuntu 9.3.0-17ubuntu1~20.04) 9.3.0, Python 3.8.5, ROOT 6.22/07
* CentOS 7, g++ (GCC) 8.3.1, Python 3.6.8, ROOT 6.20/04   
* Fedora 18, g++ (GCC) 4.7.2, Python 2.7.3, ROOT 5.34/36
* Ubuntu 12.04.5 LTS, g++ (Ubuntu/Linaro 4.6.3-1ubuntu5) 4.6.3, Python 2.7.3, ROOT 5.34/14 

## For Developers:
### Python:
```python
from specfit_cpplib import TCRFlux, TCRFluxFit, TSPECFITF1, TBPLF1, specfit_uti, specfit_canv
```
### C++:
```c++
#include "TCRFlux.h"
#include "TCRFluxFit.h"
#include "TSPECFITF1.h"
#include "TBPLF1.h"
#include "specfit_uti.h"
#include "specfit_canv.h"
```
Run
```bash
make htmldoc
```
then see 
```bash
$SPECFIT/htmldoc/index.html 
```
for additional documentation on SPECFIT functions and classes, assuming $SPECFIT is the build folder.
