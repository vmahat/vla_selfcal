import os
import subprocess
import logging
import sys
#import astropy
import matplotlib.pyplot as plt

# Define paths and parameters
msfile = sys.argv[1]  # Path to your Measurement Set
output_dir = sys.argv[2]     # Directory for outputs

"""
def findrms(mIn,maskSup=1e-7):
	
	#find the rms of an array, from Cycil Tasse/kMS

	m=mIn[np.abs(mIn)>maskSup]
	rmsold=np.std(m)
	diff=1e-1
	cut=3.
	med=np.median(m)
	for i in range(10):
		find=np.where(np.abs(m-med)<rmsold*cut)[0]
		rms=np.std(m[ind])
		if np.abs((rms-rmsold)/rmsold)<diff: 
			break
		rmsold=rms
	return rms

def flatten(f):
	#Flatten a fits file so that it becomes a 2D image. Return new header and data

	naxis=f[0].header['NAXIS']
	if naxis==2:
		return fits.PrimaryHDU(header=f[0].header,data=f[0].data)

	w = WCS(f[0].header)
	wn=WCS(naxis=2)

	wn.wcs.crpix[0]=w.wcs.crpix[0]
	wn.wcs.crpix[1]=w.wcs.crpix[1]
	wn.wcs.cdelt=w.wcs.cdelt[0:2]
	wn.wcs.crval=w.wcs.crval[0:2]
	wn.wcs.ctype[0]=w.wcs.ctype[0]
	wn.wcs.ctype[1]=w.wcs.ctype[1]

	header = wn.to_header()
	header["NAXIS"]=2
	copy=('EQUINOX','EPOCH','BMAJ', 'BMIN', 'BPA', 'RESTFRQ', 'TELESCOP', 'OBSERVER')
	for k in copy:
		r=f[0].header.get(k)
		if r is not None:
			header[k]=r

	slice=[]
	for i in range(naxis,0,-1):
		if i<=2:
			slice.append(np.s_[:],)
		else:
			slice.append(0)
		
	hdu = fits.PrimaryHDU(header=header,data=f[0].data[tuple(slice)])
	return hdu

def plotimage_aplpy(fitsimagename, outplotname, mask=None, rmsnoiseimage=None):
	import aplpy
	# image noise for plotting
	if rmsnoiseimage is None:
		hdulist = fits.open(fitsimagename)
	else:
		hdulist = fits.open(rmsnoiseimage)
		imagenoise = findrms(np.ndarray.flatten(hdulist[0].data))
		hdulist.close() 

	# image noise info
	hdulist = fits.open(fitsimagename) 
	imagenoiseinfo = findrms(np.ndarray.flatten(hdulist[0].data))
	logger.info(fitsimagename + ' Max image: ' + str(np.max(np.ndarray.flatten(hdulist[0].data))))
	logger.info(fitsimagename + ' Min image: ' + str(np.min(np.ndarray.flatten(hdulist[0].data))))
	hdulist.close()

	f = aplpy.FITSFigure(fitsimagename, slices=[0, 0])
	f.show_colorscale(vmax=16*imagenoise, vmin=-6*imagenoise, cmap='bone')
	f.set_title(fitsimagename+' (noise = {} mJy/beam)'.format(round(imagenoiseinfo*1e3, 3)))
	try: # to work around an aplpy error
		f.add_beam()
		f.beam.set_frame(True)
		f.beam.set_color('white')
		f.beam.set_edgecolor('black')
		f.beam.set_linewidth(1.)
	except:
		pass

	f.add_grid()
	f.grid.set_color('white')
	f.grid.set_alpha(0.5)
	f.grid.set_linewidth(0.2)
	f.add_colorbar()
	f.colorbar.set_axis_label_text('Flux (Jy beam$^{-1}$)')
	if mask is not None:
		try:
			f.show_contour(mask, colors='red', levels=[0.1*imagenoise], filled=False, smooth=1, alpha=0.6, linewidths=1)
		except:
			pass
	if os.path.isfile(outplotname + '.png'):
		os.system('rm -f ' + outplotname + '.png')
	f.save(outplotname, dpi=120, format='png')
	logger.info(fitsimagename + ' RMS noise: ' + str(imagenoiseinfo))
	return
"""
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
singularity_container_path = "/beegfs/general/mahatmav/lofar/long_baseline/pipeline/flocs_v5.6.0_znver2_znver2.sif" #Path to singularity container containing WSClean

casa_path = "/soft/casa-latest/bin/casa"  # Path to CASA (use the correct path)

initial_model = None               # Optional initial model

# Self-calibration parameters
#solution_intervals = ["10s","int","inf","60s","30s","inf","inf"]  # Progressive solint values
#solution_type = ["G","G","G","G","G","B","B"]
#solution_mode = ["p","p","ap","ap","ap","",""]

solution_intervals = ["inf","60s","30s","10s","int","inf","120s","inf","inf","inf"]  # Progressive solint values
solution_type = ["G","G","G","G","G","G","G","B","B","B"]
solution_mode = ["p","p","p","p","p","ap","ap","","",""]

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
		f"-no-update-model-required ",
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

# Get FITS model image from WSClean
importfits(fitsimage='{image_prefix}-MFS-model.fits', imagename='{image_prefix}-MFS-model.casaim',overwrite=True)

# Predict
ft(vis='{msfile}', model='{image_prefix}-MFS-model.casaim', usescratch=True)

# Perform gain calibration

if '{solution_mode[idx]}'=="p":
	gaincal(vis='{msfile}', caltable='{gain_table}', solint='{solint}', refant='ea23', 
	gaintype='{solution_type[idx]}', calmode='{solution_mode[idx]}')
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
