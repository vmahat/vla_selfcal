import os
import subprocess
import logging


# Set up logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger("selfcal")

# Define paths and parameters
msfile = "your_measurement_set.ms"  # Path to your Measurement Set
output_dir = "selfcal_outputs"     # Directory for outputs
os.makedirs(output_dir, exist_ok=True)

singularity_path = "/soft/singularity-3.8.4/bin/singularity" #Path to singularity installation
singularity_bind_path = "/beegfs/general/mahatmav" #Path to directory to bind
singularity_container_path = "/beegfs/general/mahatmav/lofar/long_baseline/pipeline/flocs_v5.2.0_znver2_znver2.sif" #Path to singularity container containing WSClean

casa_path = "/soft/casa-latest/bin/casa"  # Path to CASA (use the correct path)

initial_model = None               # Optional initial model

# Self-calibration parameters
solution_intervals = ["inf", "5m", "2m", "30s"]  # Progressive solint values
threshold = 0.01  # Stopping threshold for residual improvement (Jy)
gain_solutions = []  # Store gain calibration tables

# Imaging parameters for WSClean
imaging_params = {
    "size": "1024 1024",          # Image size (pixels)
    "scale": "1asec",             # Pixel scale
    "weight": "briggs 0.0",       # Weighting scheme
    "auto-threshold": "3.0",      # Threshold for CLEAN (σ)
    "auto-mask": "4.0",           # Mask threshold (σ)
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
        f"-size {imaging_params['size']}",
        f"-scale {imaging_params['scale']}",
        f"-weight {imaging_params['weight']}",
        f"-auto-threshold {imaging_params['auto-threshold']}",
        f"-auto-mask {imaging_params['auto-mask']}",
        f"-name {image_prefix}",
        msfile,
    ]
    if initial_model:
        wsclean_cmd.append(f"-model {initial_model}")

    run_command(" ".join(wsclean_cmd), shell=True)

    # CASA calibration script
    gain_table = f"{output_dir}/gains_cycle_{idx + 1}.cal"
    casa_script = f"""
from casatasks import gaincal, applycal

# Perform gain calibration
gaincal(vis='{msfile}', caltable='{gain_table}', solint='{solint}', refant='1')

# Apply calibration solutions to the MS
applycal(vis='{msfile}', gaintable={gain_solutions + [gain_table]}, calwt=False)
"""
    script_path = f"{output_dir}/casa_gaincal_cycle_{idx + 1}.py"
    with open(script_path, "w") as f:
        f.write(casa_script)

    # Run CASA script for calibration
    run_casa_command(script_path)
    gain_solutions.append(gain_table)

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
