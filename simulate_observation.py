import os
import sys
import json
import shutil
import subprocess
from astropy.coordinates import SkyCoord
import astropy.units as u
import numpy as np
from ALMA_specs import alma_dict

# === Step 1: Read config ===
if len(sys.argv) < 2:
    print("Usage: python simulate_observation.py <config.json>")
    sys.exit(1)

config_path = sys.argv[1]
with open(config_path) as f:
    config = json.load(f)

name = os.path.splitext(os.path.basename(config_path))[0]
output_dir = os.path.join("output", name)
os.makedirs(output_dir, exist_ok=True)

# === Step 2: Fill missing values from ALMA_specs ===
band = config["ALMA_band"]
config_name = config["ALMA_configuration"]
specs = alma_dict[band]

central_freq = config.get("central_freq", specs["central_freq"])
bandwidth = config.get("bandwidth", specs["bandwidth"])
pwv = config.get("pwv", specs["pwv"])
res_arcsec = specs["resolution_arcsec"][config_name]

# === Step 3: Define CASA-friendly parameters ===
project_name = f"{config['model_type']}_B{band}_{config_name}_t{config['integration_time']}s"
antennalist_file = f"alma.cycle11.{config_name.split('-')[1]}.cfg"

# === Step 4: Imaging parameters ===
skycoord = SkyCoord(config["RA"], config["Dec"], unit=(u.hourangle, u.deg))
mask_ra = skycoord.ra.to_string(unit=u.hourangle, sep='hms', precision=6, pad=True) 
mask_dec = skycoord.dec.to_string(unit=u.deg, sep='dms', precision=5, alwayssign=True)
disk_size_arcsec = config["disk_R90"] / config["distance"]  # AU
mask_semimajor = 1.5 * disk_size_arcsec 
noise_inner = 2.5 * disk_size_arcsec 
noise_outer = 3.5 * disk_size_arcsec 

cellsize = np.floor((res_arcsec / 8) * 1000) / 1000  # round down to nearest milliarcsec

if 4 * disk_size_arcsec  < 512 * cellsize:
    imsize = 512
elif 4 * disk_size_arcsec < 1024 * cellsize:
    imsize = 1024
else:
    imsize = 2048

# === Optional overrides ===
robust = config.get("robust", 0.5)
threshold = config.get("threshold", 2.0)

# === Step 5: Format CASA script ===
with open("CASA_scripts.py") as f:
    template = f.read()

filled = template.format(
    project_name=project_name,
    output_dir=output_dir,
    central_freq=central_freq,
    bandwidth=bandwidth,
    antennalist_file=antennalist_file,
    integration_time=config["integration_time"],
    pwv=pwv,
    PA=config["PA"],
    disk_R90=config["disk_R90"],
    distance=config["distance"],
    inclination=config["inclination"],
    mask_ra=mask_ra,
    mask_dec=mask_dec,
    noise_inner=noise_inner,
    noise_outer=noise_outer,
    mask_semimajor = mask_semimajor, 
    mask_semiminor = mask_semimajor * np.cos(np.radians(config["inclination"])),
    cellsize=cellsize,
    imsize=imsize,
    tclean_dir=os.path.join(output_dir, "tclean"),
    base_name=name,
    msfile=f"{project_name}/{project_name}.alma.cycle11.{config_name.split('-')[1]}.noisy.ms",
    robust=robust,
    threshold=threshold
)

# === Step 6: Save and run CASA ===
casa_script_path = os.path.join(output_dir, "simulate_tmp_casa.py")
with open(casa_script_path, "w") as f:
    f.write(filled)

subprocess.run([
    "casa", "--nogui", "--nologger", "--nologfile",
    "-c", casa_script_path
])
