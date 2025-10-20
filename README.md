#  vla_selfcal
Self-calibration pipeline to process VLA data, using CASA tools and WSClean. <br />

##  Install the prerequisite software:
-Singularity
-Singularity container from https://tikk3r.github.io/flocs/
-CASA
-Python (version 3 and above)<br />

##  Running instructions
1. Initiate a CASA environment
`module load casa` <br />

2. Ensure that PYTHONPATH points to the CASA installation python
`export PYTHONPATH="/soft/casa-latest/bin/python3"`<br />

3. Edit vla_selfcal.py so that singularity and CASA paths point to your installations. Place paths to MS and output directories <br />

3.  Run
`python3 vla_selfcal.py <ms> <output_dir> <config.txt>`<br />

## Example outputs
Basic example of running the pipeline on JVLA Ku-band C-array data on the bright radio source 3C 123, using 5 cycles with just phase calibration.
<img width="3471" height="2017" alt="vla_selfcal_3C123" src="https://github.com/user-attachments/assets/ce0679bb-ef76-47d6-9060-90ab8b252515" />
