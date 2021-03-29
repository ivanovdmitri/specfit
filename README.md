# SPECFIT

Specfit is as a Python/C++/ROOT utility for fitting arbitrary cosmic ray energy spectrum (flux) measurements into a model functions using binned Poisson log likelihood method and displaying the results.  Constant and nonlinear energy scale corrections are supported.

* Author: Dmitri Ivanov <dmiivanov@gmail.com>
* Release: 03/2021


## Easy installation and running using Dockerfile:


- Requires [Docker](https://docs.docker.com/get-docker/) installed on your system


- From the downloaded specfit project folder, run (Windows):

```shell
	C:\Users\my_user_name>\specfit docker build -t specfit .
```

OR (Linux, e.g. Ubuntu 20.04 LTS)

```shell
	my_user_name@my_host_name:~/path/to/specfit$ docker build -t specfit .
```

- Running the image will show you the manual:

```shell
	docker run --rm specfit
```

- Run the image with an interactive shell

```shell
	docker run --rm -it specfit /bin/bash
```

- Running from the inside of the container

```bash
	root@docker-desktop:/specfit/build# ./specfit.py --help
```
	    again shows the command line arguments.  Then just run

```bash
	root@docker-desktop:/specfit/build# ./specfit.py
```
	    on the command line with the desired arguments to get the desired answer.
	
- It is generally advised to run it in batch mode using option '-b'

```bash
	root@docker-desktop:/specfit/build# ./specfit.py -b
```
	    unless you have configured X11 forwarding for the Docker container.  When running in batch mode 
	    (without X11 forwarding) you can save all the plots (see specfit manual) and then copy them to 
	    some mounted volume of your choice.
		
### Show GUI on Linux: one way of configuring X11 forwarding (tested on Ubuntu 20.04 LTS) is by running the container as follows:

```bash
	$ docker run  -it --network=host --env DISPLAY=$DISPLAY  --privileged   \
	--volume="$HOME/.Xauthority:/root/.Xauthority:rw"  \
	-v /tmp/.X11-unix:/tmp/.X11-unix --rm specfit /bin/bash
```

- Then running specfit.py should produce plots in real time:

```bash
	root@docker-desktop:/specfit/build# ./specfit.py
```

### Show GUI on Windows: one way of configuring X11 forwarding (tested on Windows 10) is as follows:
- Install an X Window System Server (either one of the following would work):
  - [VcXsrv](https://sourceforge.net/projects/vcxsrv/) (free)
  - [Xming](http://www.straightrunning.com/XmingNotes/)
  - [Cygwin/X](https://x.cygwin.com/)  (free)

- Be sure that the X-server you've installed is running and accepting TCP/IP connections:
  - In the case of Xming or VcXsrv servers, choose 'Disable access control' option when running XLaunch app
  - In the case of Cygwin/X (XWin server) start it as:
  
	```bash
		$ xwin :0 -multiwindow -listen tcp
	```
- Find your IPv4 address on local network. In Windows shell, run:

```shell
	C:\Users\my_user_name>ipconfig
```
	- For example, let's assume your IPv4 address is 192.168.111.11

- Build the docker image as usual: 

```shell
	C:\Users\my_user_name>\specfit docker build -t specfit .
```
- Run the docker image as follows:

```shell
	C:\Users\my_user_name>\specfit docker run --rm -it --env DISPLAY=192.168.111.11:0.0 specfit .
```

- If everything was configured correctly, ```specfit.py`` should display plots interactively:

```bash
	root@docker-desktop:/specfit/build# ./specfit.py
```

- **Notice: in the above exercise, we are exposing the X server to the network.**  Use caution, especially if your computer is
  directly connected to the internet.

## Running on Linux without Docker

### Requirements:
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

### Install / Compile with make:
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
