# Characterisation of low-energy Compton Scattering
## Repository for the Characterisation of low-energy Compton Scattering laboratory report
All work in this repository is to be cited using the provided Citation.bib BibTeX citation

### Analysis
#### The following analysis scripts are available in the Analysis subdirectory:
* Geometry.ggb - Interavtive simulation of the experimental setup geometry used to determine uncertainty on scattering angles as well as optimum position of shielding block
* Compton Analysis.py - Analysis of the angular dependance of the energy of scattered photons
* Activity Analysis.py - Analysis of the effective activity of the source and approximate detector efficiency used for calculation of differential cross section
* Root Analysis.py - Processing of the simulation ROOT files by histogramming and writing into CSV files for further analysis
* Crossection Analysis.py - Analysis of the angular dependance of the differential cross section

### Simulation
#### The following steps are required to run the GEANT4 simulation:
1. change the CMake prefix path to the path of your GEANT4 installation in the `Simulation/source/CmakeLists.txt` file
   ```cmake:
   list(APPEND CMAKE_PREFIX_PATH "/usr/local/Geant4/geant4-v11.1.3-install/lib/cmake/Geant4")
   ```
3. create a build directory and cd into it
4. run CMake with the following command
   ```zsh:
   cmake ../source
   ```
6. build the simulation with the make command where `-j` specifies the number of thread to use
   ```zsh:
   make -j4
   ```
8. run the simulation in visualisation mode with the preset geometry
   ```zsh:
   ./scat
   ```
10. now run the simulation in batch mode using the `scat.mac` macro file
   ```zsh:
   ./scar scat.mac
   ```

#### The following messenger commands can be used to control the simulation
* `/setup/targ <bool>` set to construct construct the scattering target
* `/setup/shield <bool>` set to construct the shielding block
* `/setup/block <bool>` set to block the source for background measurements
* `/setup/scatAng <double>` specify the scattering angle at which to position the source
* `/setup/shieldAng <double>` specify the shielding angle at which to position the shielding block
* `/setup/shieldDist <double>` specify the distance from the center of the target to the front face of the shielding block
  
* `/setup/partE <double>` specify energy of simulated photons in keV
* `/setup/momRand <bool>` set to generate photons with random momentum direction vector
* `/setup/momRest <bool>` set to restrict the momentum direction vector to increase simulation statistic has to be accounted for by factor in Root Analysis.py
* `/setup/momX <double>` set x component of momentum direction vector before source is rotated to specified scattering angle can also set `momY` and `momZ`

### Dependencies
#### The following dependencies are required to run the complete simulation and analysis:
* [GEANT4](https://geant4.web.cern.ch)
* [Qt5](https://doc.qt.io/qt.html#qt5)
* [CMake](https://cmake.org)
* [UpRoot](https://uproot.readthedocs.io/en/latest/)

Lukas Kostal, 24.11.2023, ICL
