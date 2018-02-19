# boost-geometry-proposal
This repository contains the programming competency test and GSoC proposal.

## GeoLib

### Build instructions
First clone the repository with the following command:
```bash
$ git clone https://github.com/adl1995/boost-geometry-proposal.git
```
Then, `cd` into the cloned repository by:
```bash
$ cd boost-geometry-proposal
```
To compile the library, first create a build directory:
```bash
$ mkdir build
$ cd build
```
The next step is to run CMake to configure the project:
```bash
$ cmake ../
```
Once CMake is configured, the library can be built by typing `make`. This will build the 'geolib_tests' component:
```bash
$ make
```
The tests can be run by typing:
```bash
$ make tests
```
To build the documentation using Doxygen, type:
```bash
$ make docs
```
Finally, the HTML documentation can be opened with Firefox by typing:
```bash
$ firefox docs/html/index.html
```
