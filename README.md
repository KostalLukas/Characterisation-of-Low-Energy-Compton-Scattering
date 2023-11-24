# Characterisation of low-energy Compton Scattering
## Repository for the Characterisation of low-energy Compton Scattering laboratory report
All work in this repository is to be cited using the provided Citation.bib BibTeX citation

### Analysis
The following analysis scripts are available in the Analysis subdirectory:
* Geometry.ggb - Interavtive simulation of the experimental setup geometry used to determine uncertainty on scattering angles as well as optimum position of shielding block
* Compton Analysis.py - Analysis of the angular dependance of the energy of scattered photons
* Activity Analysis.py - Analysis of the effective activity of the source and approximate detector efficiency used for calculation of differential cross section
* Root Analysis.py - Processing of the simulation ROOT files by histogramming and writing into CSV files for further analysis
* Crossection Analysis.py - Analysis of the angular dependance of the differential cross section

### Simulation
The following steps are required to run the GEANT4 simulation:

### Dependancies
The following dependancies are required to run the complete simulation and analysis:
* [GEANT4](https://geant4.web.cern.ch)
* [Qt5](https://doc.qt.io/qt.html#qt5)
* [CMake](https://cmake.org)
* [UpRoot](https://uproot.readthedocs.io/en/latest/)

Lukas Kostal, 24.11.2023, ICL
