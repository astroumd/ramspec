
Below are steps to make spectrum curves from a RAMSES output. This is tested on my star formation simulations. 

1. Download ramspec (this repo) including its submodules: `git clone --recursive https://github.com/astroumd/ramspec.git`
3. Compile ramski:
   1. `cd ramspec/external/RAMSKI_test/f90; make`
   2. Note: check if `gfortran` exists via `which gfortran`. You may need to load the necessary module, e.g. `module load gcc` (on HPC) or `brew install gcc` (on macOS). 
4. Compile amr2cube: `cd ramspec/tools/amr2cube; make`
5. Install SKIRT following the official website. Install PTS as well and add it to PYTHONPATH (see tests/run.sh for how to do this). Confirm by running `python -c 'import pts'`
6. Install necessary packages for pts
   1. List Python package dependencies: `python -m pts.do list_dependencies` 
   2. Install missing packages via conda or pip. For example, I have to do the following: `conda install lxml reportlab` and `python -m pip install fsps opencv-python`
7. Prepare a namelist file for RAMSKI/SKIRT simulation and run the test:
   1. Copy the tests folder to somewhere in your computer. 
   2. Set PYTHONPATH in your shell following the example in run.sh. You should set your PYTHONPATH to the directory containing ramspec AND the path to PTS.
   3. Run import_modules.py to make sure all modules and installed and correctly linked: `cd tests; python test_modules.py`
   4. Run the test: `python test_run.py`. Before running this test, make sure to change a few variables in test_run.py: `ramspec.FIELDS`, `jobpath`, `nml`, `outs`. 

Notes:
1. Use ramspec as a Python module and don't change its contents unless you know what you are doing. Your change will conflict with future git pull as I keep updating this repository. Make a folder outside of ramspec, copy test_run.py and do your run there. 
