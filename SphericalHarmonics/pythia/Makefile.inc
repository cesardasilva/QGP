# PYTHIA configuration file.
# Generated on Mon Nov 7 16:28:05 EST 2022 with the user supplied options:
# --prefix=/cvmfs/sphenix.sdcc.bnl.gov/gcc-12.1.0/opt/sphenix/core/pythia8307_hepmc3-3.2.5_root-6.26.06
# --with-hepmc3=/cvmfs/sphenix.sdcc.bnl.gov/gcc-12.1.0/opt/sphenix/core/HepMC3-3.2.5_root-6.26.06
# --with-hepmc3-lib=/cvmfs/sphenix.sdcc.bnl.gov/gcc-12.1.0/opt/sphenix/core/HepMC3-3.2.5_root-6.26.06/lib64
# --with-evtgen=/cvmfs/sphenix.sdcc.bnl.gov/gcc-12.1.0/opt/sphenix/core/EvtGen-R02-02-00_hepmc3-3.2.5_root-6.26.06
# --with-evtgen-lib=/cvmfs/sphenix.sdcc.bnl.gov/gcc-12.1.0/opt/sphenix/core/EvtGen-R02-02-00_hepmc3-3.2.5_root-6.26.06/lib64
# --with-hepmc2=/cvmfs/sphenix.sdcc.bnl.gov/gcc-12.1.0/opt/sphenix/core/HepMC-2.06.11
# --with-lhapdf6=/cvmfs/sphenix.sdcc.bnl.gov/gcc-12.1.0/opt/sphenix/core/LHAPDF-6.5.3
# --with-fastjet3=/cvmfs/sphenix.sdcc.bnl.gov/gcc-12.1.0/opt/sphenix/core/fastjet-3.4.0
# --with-python=/cvmfs/sphenix.sdcc.bnl.gov/gcc-12.1.0/opt/sphenix/core/include/python3.8
# --enable-optdebug
# --with-gzip
# --with-python=/cvmfs/sphenix.sdcc.bnl.gov/gcc-12.1.0/opt/sphenix/core/include/python3.10

# Install directory prefixes.
PREFIX_BIN=/cvmfs/sphenix.sdcc.bnl.gov/gcc-12.1.0/opt/sphenix/core/pythia8307_hepmc3-3.2.5_root-6.26.06/bin
PREFIX_INCLUDE=/cvmfs/sphenix.sdcc.bnl.gov/gcc-12.1.0/opt/sphenix/core/pythia8307_hepmc3-3.2.5_root-6.26.06/include
PREFIX_LIB=/cvmfs/sphenix.sdcc.bnl.gov/gcc-12.1.0/opt/sphenix/core/pythia8307_hepmc3-3.2.5_root-6.26.06/lib
PREFIX_SHARE=/cvmfs/sphenix.sdcc.bnl.gov/gcc-12.1.0/opt/sphenix/core/pythia8307_hepmc3-3.2.5_root-6.26.06/share/Pythia8

# Compilation flags (see ./configure --help for further documentation).
CXX=/cvmfs/sphenix.sdcc.bnl.gov/gcc-12.1.0/opt/sphenix/core/gcc/12.1.0-57c96/x86_64-centos7/bin/g++
CXX_COMMON=-g -O -std=c++11 -pedantic -W -Wall -Wshadow -fPIC -pthread -DGZIP
CXX_SHARED=-shared
CXX_SONAME=-Wl,-soname,
LIB_SUFFIX=.so
OBJ_COMMON=

EVTGEN_USE=true
EVTGEN_CONFIG=
EVTGEN_BIN=/cvmfs/sphenix.sdcc.bnl.gov/gcc-12.1.0/opt/sphenix/core/EvtGen-R02-02-00_hepmc3-3.2.5_root-6.24.06/
EVTGEN_INCLUDE=-I/cvmfs/sphenix.sdcc.bnl.gov/gcc-12.1.0/opt/sphenix/core/EvtGen-R02-02-00_hepmc3-3.2.5_root-6.24.06/include -I/cvmfs/sphenix.sdcc.bnl.gov/gcc-12.1.0/opt/sphenix/core/HepMC-2.06.11/include
EVTGEN_LIB=-L/cvmfs/sphenix.sdcc.bnl.gov/gcc-12.1.0/opt/sphenix/core/EvtGen-R02-02-00_hepmc3-3.2.5_root-6.24.06/lib64 -Wl,-rpath,/cvmfs/sphenix.sdcc.bnl.gov/gcc-12.1.0/opt/sphenix/core/EvtGen-R02-02-00_hepmc3-3.2.5_root-6.24.06/lib64 -lEvtGen -lEvtGenExternal -L/cvmfs/sphenix.sdcc.bnl.gov/gcc-12.1.0/opt/sphenix/core/HepMC-2.06.11/lib -Wl,-rpath,/cvmfs/sphenix.sdcc.bnl.gov/gcc-12.1.0/opt/sphenix/core/HepMC-2.06.11/lib -lHepMC

FASTJET3_USE=true
FASTJET3_CONFIG=fastjet-config
FASTJET3_BIN=/cvmfs/sphenix.sdcc.bnl.gov/gcc-12.1.0/opt/sphenix/core/fastjet-3.4.0/bin/
FASTJET3_INCLUDE=-I/cvmfs/sphenix.sdcc.bnl.gov/gcc-12.1.0/opt/sphenix/core/fastjet-3.4.0/include
FASTJET3_LIB=-L/cvmfs/sphenix.sdcc.bnl.gov/gcc-12.1.0/opt/sphenix/core/fastjet-3.4.0/lib -Wl,-rpath,/cvmfs/sphenix.sdcc.bnl.gov/gcc-12.1.0/opt/sphenix/core/fastjet-3.4.0/lib -lfastjet

HEPMC2_USE=true
HEPMC2_CONFIG=
HEPMC2_BIN=/cvmfs/sphenix.sdcc.bnl.gov/gcc-12.1.0/opt/sphenix/core/HepMC-2.06.11/
HEPMC2_INCLUDE=-I/cvmfs/sphenix.sdcc.bnl.gov/gcc-12.1.0/opt/sphenix/core/HepMC-2.06.11/include
HEPMC2_LIB=-L/cvmfs/sphenix.sdcc.bnl.gov/gcc-12.1.0/opt/sphenix/core/HepMC-2.06.11/lib -Wl,-rpath,/cvmfs/sphenix.sdcc.bnl.gov/gcc-12.1.0/opt/sphenix/core/HepMC-2.06.11/lib -lHepMC

HEPMC3_USE=true
HEPMC3_CONFIG=HepMC3-config
HEPMC3_BIN=/cvmfs/sphenix.sdcc.bnl.gov/gcc-12.1.0/opt/sphenix/core/HepMC3-3.2.5_root-6.26.06/bin/
HEPMC3_INCLUDE=-I/cvmfs/sphenix.sdcc.bnl.gov/gcc-12.1.0/opt/sphenix/core/HepMC3-3.2.5_root-6.26.06/include
HEPMC3_LIB=-L/cvmfs/sphenix.sdcc.bnl.gov/gcc-12.1.0/opt/sphenix/core/HepMC3-3.2.5_root-6.26.06/lib64 -Wl,-rpath,/cvmfs/sphenix.sdcc.bnl.gov/gcc-12.1.0/opt/sphenix/core/HepMC3-3.2.5_root-6.26.06/lib64 -lHepMC3

LHAPDF5_USE=false
LHAPDF5_CONFIG=
LHAPDF5_BIN=
LHAPDF5_INCLUDE=
LHAPDF5_LIB=

LHAPDF6_USE=true
LHAPDF6_CONFIG=lhapdf-config
LHAPDF6_BIN=/cvmfs/sphenix.sdcc.bnl.gov/gcc-12.1.0/opt/sphenix/core/LHAPDF-6.5.3/bin/
LHAPDF6_INCLUDE=-I/cvmfs/sphenix.sdcc.bnl.gov/gcc-12.1.0/opt/sphenix/core/LHAPDF-6.5.3/include
LHAPDF6_LIB=-L/cvmfs/sphenix.sdcc.bnl.gov/gcc-12.1.0/opt/sphenix/core/LHAPDF-6.5.3/lib -Wl,-rpath,/cvmfs/sphenix.sdcc.bnl.gov/gcc-12.1.0/opt/sphenix/core/LHAPDF-6.5.3/lib -lLHAPDF

POWHEG_USE=false
POWHEG_CONFIG=
POWHEG_BIN=
POWHEG_INCLUDE=
POWHEG_LIB=

RIVET_USE=false
RIVET_CONFIG=
RIVET_BIN=
RIVET_INCLUDE=
RIVET_LIB=

ROOT_USE=true
ROOT_CONFIG=root-config
ROOT_BIN=/cvmfs/sphenix.sdcc.bnl.gov/gcc-12.1.0/opt/sphenix/core/root-6.26.06.p01/bin/
ROOT_INCLUDE=/cvmfs/sphenix.sdcc.bnl.gov/gcc-12.1.0/opt/sphenix/core/root-6.26.06.p01/include
ROOT_LIB=-L/cvmfs/sphenix.sdcc.bnl.gov/gcc-12.1.0/opt/sphenix/core/root-6.26.06.p01/lib -lCore -lImt -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lROOTVecOps -lTree -lTreePlayer -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread -lMultiProc -lROOTDataFrame -Wl,-rpath,/cvmfs/sphenix.sdcc.bnl.gov/gcc-12.1.0/opt/sphenix/core/root-6.26.06.p01/lib -pthread -lm -ldl -rdynamic

GZIP_USE=true
GZIP_CONFIG=
GZIP_BIN=
GZIP_INCLUDE=
GZIP_LIB=-lz

PYTHON_USE=true
PYTHON_CONFIG=python-config
PYTHON_BIN=/cvmfs/sphenix.sdcc.bnl.gov/gcc-12.1.0/opt/sphenix/core/stow/Python-3.10.8/include/python3.10/
PYTHON_INCLUDE=-I/cvmfs/sphenix.sdcc.bnl.gov/gcc-12.1.0/opt/sphenix/core/stow/Python-3.10.8/include/python3.10
PYTHON_LIB=-L/cvmfs/sphenix.sdcc.bnl.gov/gcc-12.1.0/opt/sphenix/core/stow/Python-3.10.8/include/python3.10 -Wl,-rpath,/cvmfs/sphenix.sdcc.bnl.gov/gcc-12.1.0/opt/sphenix/core/stow/Python-3.10.8/include/python3.10

MG5MES_USE=false
MG5MES_CONFIG=
MG5MES_BIN=
MG5MES_INCLUDE=
MG5MES_LIB=

OPENMP_USE=false
OPENMP_CONFIG=
OPENMP_BIN=
OPENMP_INCLUDE=
OPENMP_LIB=

MPICH_USE=false
MPICH_CONFIG=
MPICH_BIN=
MPICH_INCLUDE=
MPICH_LIB=

HDF5_USE=false
HDF5_CONFIG=
HDF5_BIN=
HDF5_INCLUDE=
HDF5_LIB=

HIGHFIVE_USE=false
HIGHFIVE_CONFIG=
HIGHFIVE_BIN=
HIGHFIVE_INCLUDE=
HIGHFIVE_LIB=
