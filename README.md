# UMLINCA (LinThurber)

[![License](https://img.shields.io/badge/License-BSD_3--Clause-blue.svg)](https://opensource.org/licenses/BSD-3-Clause)
![GitHub repo size](https://img.shields.io/github/repo-size/sceccode/uwlinca)
[![uwlinca-ucvm-ci Actions Status](https://github.com/SCECcode/uwlinca/workflows/uwlinca-ucvm-ci/badge.svg)](https://github.com/SCECcode/uwlinca/actions)

The UW statewide model is a three-dimensional (3D) tomographic model of the P wave velocity (Vp) structure 
of northern California. It was obtained using a regional-scale double-difference tomography algorithm that 
incorporates a finite-difference travel time calculator and spatial smoothing constraints. Arrival times 
from earthquakes and travel times from controlled-source explosions, recorded at network and/or temporary
stations, were inverted for Vp on a 3D grid with horizontal node spacing of 10 to 20 km and vertical node 
spacing of 3 to 8. 

## Installation

This package is intended to be installed as part of the UCVM framework,
version 25.x or higher.

## Contact the authors

If you would like to contact the authors regarding this software,
please e-mail software@scec.org. Note this e-mail address should
be used for questions regarding the software itself (e.g. how
do I link the library properly?). Questions regarding the model's
science (e.g. on what paper is the UWLINCA based?) should be directed
to the model's authors, located in the AUTHORS file.

## NOTE

A right rectangle with 36 degree rotation

vs and density are calculated, from https://pubs.usgs.gov/of/2005/1317/of2005-1317.pdf

<pre>
  *[eqn. 1] Vs (km/s) = 0.7858 – 1.2344Vp + 0.7949Vp2 – 0.1238Vp3 + 0.0064Vp4
  *[eqn. 6] r (g/cm3) = 1.6612Vp – 0.4721Vp2 + 0.0671Vp3 – 0.0043Vp4 + 0.000106Vp5

</pre>

