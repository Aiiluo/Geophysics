![中国石油大学（华东）](logo_UPC.jpg)


# Geophysical-data-processing-methods
Geophysical data processing methods, Program code implementation of geophysical data processing methods, including parallel algorithms (cuda, mpich, openmp, etc.) with finite difference forward modeling, inverse time migration, and full waveform inversion.
* UPC: <http://www.upc.edu.cn/>
* my blog: <https://blog.csdn.net/Rong_Toa>

# Key Words 

## Geophysical Medium
* 2D (two dimension)
* 3D (three dimension)
* VTI (a kind of anisotropic）
* TTI (a kind of anisotropic)
* ISO (isotropic)

## Geophysical Wave
* acoustic wave (P, qP, "qSV")
* Elastic wave (P, SH, SV)

## Geophysical Methods
* FD (Finite Difference)
* RTM (Reverse Time Migration)
* FWI (Full Waveform Inversion)
* TOMO (Tomography)
* RayTracing

## Geophysical Gather
* ADCIGs (Angel Domain Common Imaging Gather)

## Geophysical file format
* binary/data 
* SU (Seismic Unix.su)
* segy (seg-Y .segy or .sgy)

## Tools & Software
* HPC (High-performance computing)
* omp (openmp, multi-thread of CPU)
* mpi (mpich or openmpi, multi-core/node)
* cuda (NVIDIA CUDA-Toolkit, version>=7.5, multi-thread of GPU)

## GUI Tools
* GUI (Graphical User Interface)
* gtk (GIMP Toolkit, version>=2.0, such as "gtk+-2.0" or "gtk+-3.0")
* Qt (Qt-Creater Toolkit, C++)
* software

## Others
* utils （tools function)
* mix (mix programming, mix-mpi+cuda meas mix mpi and cuda)
* CWP (Seismic Unix:<http://www.seismicunix.com/w/Main_Page>)
* Madagascar (GitHub src Link:<https://github.com/ahay/src>)

## Detail in Each README.md in Each Directory

## Dependence && Envrionment
### OS
* linux
### Compiler
* gcc, nvcc, mpicc, javac
### Software(some)
* gcc, cuda, mpich/openmpi, openmp, Qt-Creater, gtk+-2.0/3.0, JDK, make/cmake
## GUI Software in This Project
### Software 1:
* in directory "software-2D-ISO-VTI-FD-RayTracing-javaSwingAwt-GUI-0.0" is a software programming with java&c, you can see the detail in the README file in this directory.
### Software 2:
* in dirctory "2D-3D-vti-fd-cuda-gtk-GUI" is a software programming with c+cuda+gtk, you can check detail in README file in that directory.

