import os
import subprocess
import logging
import sys
#import astropy
import matplotlib.pyplot as plt
import ast
import csv
import shutil

def load_calibration_config(config_path):
    config = {}
    with open(config_path, "r") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue  # Skip empty lines and comments
            if "=" not in line:
                continue
            key, value = line.split("=", 1)
            key = key.strip()
            value = value.strip()
            try:
                config[key] = ast.literal_eval(value)
            except Exception as e:
                raise ValueError(f"Failed to parse line: {line}\n{e}")
    return config

if len(sys.argv) < 4:
    print("Usage: python selfcal.py <msfile> <output_dir> <config_file.py>")
    sys.exit(1)

msfile = sys.argv[1]
output_dir = sys.argv[2]
config_file = sys.argv[3]

config = load_calibration_config(config_file)

# Extract the lists
solution_intervals = config.get("solution_intervals")
solution_type = config.get("solution_type")
solution_mode = config.get("solution_mode")
flag_residual = config.get("flag_residual")
smooth_solutions = config.get("smooth_solutions")

# Optional: Validate lengths
if not (len(solution_intervals) == len(solution_type) == len(solution_mode)):
    print("Mismatch in lengths of solution lists.")
    sys.exit(1)

# Set up logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger("selfcal")


try:
	os.makedirs(output_dir, exist_ok=True)
except:
	print("Cannot make output directory!")
	sys.exit()

singularity_path = "/soft/singularity-3.8.4/bin/singularity" #Path to singularity installation
singularity_bind_path = "/beegfs/general/mahatmav" #Path to directory to bind
singularity_container_path = "/beegfs/general/mahatmav/lofar/long_baseline/pipeline/flocs_v5.7.1_znver2_znver2.sif" #Path to singularity container containing WSClean

casa_path = "/soft/casa-latest/bin/casa"  # Path to CASA (use the correct path)

initial_model = None               # Optional initial model

with open(config_file, 'r') as f:
    reader = csv.reader(f)
    for row in reader:
        if not row or row[0].startswith("#"):
            continue  # skip comments/empty lines
        solint = row[0].strip()
        soltype = row[1].strip() if len(row) > 1 else ""
        solmode = row[2].strip() if len(row) > 2 else ""
        solution_intervals.append(solint)
        solution_type.append(soltype)
        solution_mode.append(solmode)

threshold = 0.01  # Stopping threshold for residual improvement (Jy)
gain_solutions = []  # Store gain calibration tables
continue_imaging=False #Option to re-run imaging even if images exist
# Imaging parameters for WSClean
imaging_params = {
	"size": "4096 4096",          # Image size (pixels)
	"scale": "0.02asec",             # Pixel scale
	"weight": "briggs -0.5",       # Weighting scheme
	"auto-threshold": "0.5",      # Threshold for CLEAN (σ)
	"auto-mask": "3.0",           # Mask threshold (σ)
	"niter": "500000",			# Minor cycles
	"nmiter": "50",
	"mgain": "0.7",				# Major cycles
	"channels-out": "12",		# Channels to do peak-finding
	"fit-spectral-pol": "4",		# MFS	
	"padding": "1.4",
	"max-scales": "8"
}

# Function to run a command and check output
def run_command(cmd, shell=False):
	try:
		logger.info(f"Running: {' '.join(cmd) if isinstance(cmd, list) else cmd}")
		subprocess.run(cmd, shell=shell, check=True)
	except subprocess.CalledProcessError as e:
		logger.error(f"Command failed: {e}")
		raise

# Function to run CASA commands
def run_casa_command(casa_script):
	casa_cmd = f"{casa_path} --nogui -c {casa_script}"
	run_command(casa_cmd, shell=True)

#Make a copy of the MS in the output dir
if os.path.exists(os.path.join(output_dir,os.path.basename(msfile))):
	pass
else:
	shutil.copytree(msfile,os.path.join(output_dir,os.path.basename(msfile)))
	msfile=os.path.join(output_dir,os.path.basename(msfile))

# Self-calibration loop
for idx, solint in enumerate(solution_intervals):

	logger.info(f"Starting self-calibration cycle {idx + 1} with solint={solint}")

    # WSClean imaging step
	image_prefix = f"{output_dir}/selfcal_cycle_{idx + 1}"

	wsclean_cmd = [
		singularity_path,
		"run",
		f"--bind {singularity_bind_path} {singularity_container_path}",
		"wsclean",
		f"-update-model-required ",
		f"-reorder ",
		f"-size {imaging_params['size']}",
		f"-scale {imaging_params['scale']}",
		f"-weight {imaging_params['weight']}",
		f"-auto-threshold {imaging_params['auto-threshold']}",
		f"-auto-mask {imaging_params['auto-mask']}",
		f"-name {image_prefix}",
		f"-niter {imaging_params['niter']}",
		f"-nmiter {imaging_params['nmiter']}",
		f"-mgain {imaging_params['mgain']}",
		f"-channels-out {imaging_params['channels-out']}",
		f"-fit-spectral-pol {imaging_params['fit-spectral-pol']}",
		f"-join-channels ",
		f"-multiscale ",
		f"-multiscale-max-scales {imaging_params['max-scales']}",
		f"-multiscale-scale-bias 0.8",
		f"-padding {imaging_params['padding']}",
		f"-local-rms ",
		f"-local-rms-strength 0.5",
		msfile,
	]
	#if initial_model: This doesn't need to exist, surely
	#    wsclean_cmd.append(f"-model {initial_model}")
	if continue_imaging and os.path.exists(f"{output_dir}/selfcal_cycle_{idx + 1}-MFS-image.fits"): #don't image if exists
		print(f'First image exists already! Proceeding with calibration')
	else:
		run_command(" ".join(wsclean_cmd), shell=True)
		"""
		if {imaging_params['channels-out']}>1:
			fitsimage={image_prefix}+"-MFS-image.fits"
		pngimage=f"{output_dir}/selfcal_cycle_{idx + 1}"
		plotimage_aplpy(fitsimage,pngimage)
		"""
	# CASA calibration script
	gain_table = f"{output_dir}/gains_cycle_{idx + 1}.cal"
	casa_script = f"""
from casatasks import gaincal, applycal

#store some global variables in this local script
solution_mode = {solution_mode}
solution_intervals = {solution_intervals}
solution_type = {solution_type}
temp_gain_table = None
temp_p_gain_table=None
temp_ap_gain_table=None

# Perform gain calibration

if '{solution_mode[idx]}'=="p":
	gaincal(vis='{msfile}', caltable='{gain_table}', solint='{solint}', refant='ea23', 
	gaintype='{solution_type[idx]}', calmode='{solution_mode[idx]}')
	if "{smooth_solutions}":
		smoothcal(vis="{msfile}",tablein='{gain_table}',smoothtype='mean',
		smoothtime=10.0)
if '{solution_mode[idx]}'=="ap":
	#Check the last phase-only and find its solint to do another round to pre-apply to ap
	for prev_idx in range({idx}-1,-1,-1): #loop backwards from previous to first
		if solution_mode[prev_idx] == "p":
			prev_solint=solution_intervals[prev_idx]
			break
	temp_gain_table = f"{output_dir}/temp_pre_ap_cycle{idx+1}.cal"
	gaincal(vis='{msfile}', caltable=temp_gain_table, solint=prev_solint, refant='ea23', 
		gaintype=solution_type[prev_idx], calmode=solution_mode[prev_idx])
	#Now do ap with previous p on the fly
	gaincal(vis='{msfile}', caltable='{gain_table}', solint='{solint}', refant='ea23', 
		gaintype='{solution_type[idx]}', calmode='{solution_mode[idx]}', 
		gaintable=[temp_gain_table],solnorm=True)
	if "{smooth_solutions}":
		smoothcal(vis="{msfile}",tablein='{gain_table}',smoothtype='mean',
		smoothtime=60.0)
if '{solution_type[idx]}'=="B":
		#Check the last phase-only and find its solint to do another round to pre-apply to ap
	for prev_idx in range({idx}-1,-1,-1): #loop backwards from previous to first
		if solution_mode[prev_idx] == "p":
			prev_solint=solution_intervals[prev_idx]
			break
	temp_p_gain_table = f"{output_dir}/temp_pre_bp_p_cycle{idx+1}.cal"
	gaincal(vis='{msfile}', caltable=temp_p_gain_table, solint=prev_solint, refant='ea23', 
		gaintype=solution_type[prev_idx], calmode=solution_mode[prev_idx])
	#Now do ap with previous p on the fly
	for prev_idx in range({idx}-1,-1,-1): #loop backwards from previous to first
		if solution_mode[prev_idx] == "ap":
			prev_solint=solution_intervals[prev_idx]
			break
	temp_ap_gain_table = f"{output_dir}/temp_pre_bp_ap_cycle{idx+1}.cal"
	gaincal(vis='{msfile}', caltable=temp_ap_gain_table, solint=prev_solint, refant='ea23', 
		gaintype=solution_type[prev_idx], calmode=solution_mode[prev_idx], 
		gaintable=[temp_p_gain_table],solnorm=True)

	#Now do bp with previous p and ap solutions
	bandpass(vis='{msfile}', caltable='{gain_table}', solint='{solint}', refant='ea23', 
		bandtype='{solution_type[idx]}', 
		gaintable=[temp_p_gain_table,temp_ap_gain_table],solnorm=False)
	if "{smooth_solutions}":
		smoothcal(vis="{msfile}",tablein='{gain_table}',smoothtype='mean',
		smoothtime=600.0)

#Flag solutions
if '{flag_residual}':
	print('Using rflag algorithm to flag outlier solutions...')
	flagdata(vis='{gain_table}',mode='rflag',datacolumn='CPARM')

# Collect recently created tables

current_gain_tables = []
if temp_gain_table:
	current_gain_tables.append(temp_gain_table)#should include latest phase only
if temp_p_gain_table:
	current_gain_tables.append(temp_p_gain_table)#should include latest phase only
if temp_ap_gain_table:
	current_gain_tables.append(temp_ap_gain_table)#should include latest a+phase only
current_gain_tables.append('{gain_table}')

# Apply calibration solutions to the MS
applycal(vis='{msfile}', gaintable=current_gain_tables, calwt=False)
	"""
	script_path = f"{output_dir}/casa_gaincal_cycle_{idx + 1}.py"
	with open(script_path, "w") as f:
		f.write(casa_script)

    # Run CASA script for calibration
	run_casa_command(script_path)
	#gain_solutions.append(gain_table) Don't apply previous solutions (only needed when doing amp selfcal)

	# Check if improvement is sufficient
	residual_image = f"{image_prefix}-residual.fits"
	if idx > 0:
		prev_residual_image = f"{output_dir}/selfcal_cycle_{idx}-residual.fits"
		# Compare residuals here (use a tool to assess noise or improvement threshold)
		logger.info(f"Compare {residual_image} with {prev_residual_image}")
		# If no significant improvement, break the loop

    # Prepare for next cycle
	initial_model = f"{image_prefix}-model.fits"

logger.info("Self-calibration completed successfully!")
