#  vla_selfcal
Self-calibration pipeline to process VLA data, using CASA tools and WSClean.

##  Install the prerequisite software:
-Singularity
-Singularity container from https://tikk3r.github.io/flocs/
-CASA
-Python (version 3 and above)

##  Running instructions
Initiate a CASA environment
`module load casa`

### Ensure that PYTHONPATH points to the CASA installation python
`export PYTHONPATH="/soft/casa-latest/bin/python3"`

### Edit vla_selfcal.py so that singularity and CASA paths point to your installations. Place paths to MS and output directories

###  Run
`python3 vla_selfcal.py`
