
1. Download ramspectrum (this repo) including its submodules:
   1. `git clone --recursive https://github.com/astroumd/ramspectrum.git`
2. Compile ramski:
   1. `cd ramspectrum/external/RAMSKI_test/f90; make`
   2. Note: check if `gfortran` exists via `which gfortran`. You may need to load the necessary module, e.g. `module load gcc`
3. Compile amr2cube:
   1. `cd ramspectrum/tools/amr2cube; make`
4. Install SKIRT following the official website. Install PTS as well and add it to PYTHONPATH. Check it by running `python -c 'import pts'`
5. Check and install necessary packages for pts
   1. List Python package dependencies: `python -m pts.do list_dependencies` and install missing packages via conda or pip. For example, I have to do the following: `conda install lxml reportlab` and `python -m pip install fsps opencv-python`
6. Run the tests
   1. First, run the import_modules test to make sure all modules and installed and correctly linked: `cd test; python test_modules.py`
   2. Run the actual test: `python test_run.py`. Before running this test, make sure to change a few variables: `ramspec.FIELDS`, `jobpath`, `nml`, `outs`, 

Assumptions:
1. The sink particle data is the same as my star formation runs. This is probably not always true. 


