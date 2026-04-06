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
bandwidth    = config.get("bandwidth",    specs["bandwidth"])
pwv          = config.get("pwv",          specs["pwv"])
res_arcsec   = specs["resolution_arcsec"][config_name]

# === Step 3: CASA-friendly parameters ===
project_name    = f"{config['model_type']}_B{band}_{config_name}_t{config['integration_time']}s"
antennalist_file = f"alma.cycle11.{config_name.split('-')[1]}.cfg"

# === Step 4: Imaging parameters ===
skycoord       = SkyCoord(config["RA"], config["Dec"], unit=(u.hourangle, u.deg))
mask_ra        = skycoord.ra.to_string(unit=u.hourangle, sep='hms', precision=6, pad=True)
mask_dec       = skycoord.dec.to_string(unit=u.deg, sep='dms', precision=5, alwayssign=True)
disk_size_arcsec = config["disk_R90"] / config["distance"]

mask_semimajor = 1.5 * disk_size_arcsec
mask_semiminor = mask_semimajor * np.cos(np.radians(config["inclination"]))
noise_inner    = 2.5 * disk_size_arcsec
noise_outer    = 3.5 * disk_size_arcsec

# Check for accidental underscore-prefixed override keys
for k in ["_cellsize", "_imsize"]:
    if k in config:
        print(f"  NOTE: found \"{k}\" in config — this is INACTIVE. "
              f"Rename to \"{k[1:]}\" to activate the manual override.")

if "cellsize" in config:
    cellsize = config["cellsize"]
    print(f"  cellsize: MANUAL OVERRIDE = {cellsize:.4f} arcsec")
else:
    cellsize = max(np.floor((res_arcsec / 8) * 1000) / 1000, 0.002)
    print(f"  cellsize: AUTO = {cellsize:.4f} arcsec  (beam/8 = {res_arcsec/8:.4f})")

if "imsize" in config:
    imsize = config["imsize"]
    print(f"  imsize: MANUAL OVERRIDE = {imsize}")
else:
    fov_needed = 2 * noise_outer * 1.2
    if fov_needed < 512 * cellsize:
        imsize = 512
    elif fov_needed < 1024 * cellsize:
        imsize = 1024
    elif fov_needed < 2048 * cellsize:
        imsize = 2048
    elif fov_needed < 4096 * cellsize:
        imsize = 4096
    else:
        imsize = 8192
        print(f"  WARNING: very large image ({imsize} px). "
              f"Consider setting \"imsize\" and \"cellsize\" manually in the config.")
    print(f"  imsize: AUTO = {imsize}  (FOV needed: {fov_needed:.3f}\"  covers: {imsize*cellsize:.3f}\")")


# ------------------------------------------------------------------ #
# Geometry validation: noise annulus must fit inside the image
# ------------------------------------------------------------------ #
image_half_width = (imsize * cellsize) / 2.0   # arcsec from centre to edge
if noise_outer >= image_half_width:
    print()
    print("  *** GEOMETRY ERROR ***")
    print(f"  Noise annulus outer radius : {noise_outer:.3f} arcsec")
    print(f"  Image half-width           : {image_half_width:.3f} arcsec  ({imsize} px x {cellsize:.4f} arcsec/px)")
    print(f"  The noise annulus falls OUTSIDE the image.")
    print(f"  Fix: increase imsize, increase cellsize, or reduce noise_outer.")
    print(f"  Suggestion: imsize >= {int(2 * noise_outer / cellsize) + 1} px at current cellsize.")
    import sys; sys.exit(1)
elif noise_outer >= image_half_width * 0.9:
    print(f'  WARNING: noise annulus outer edge ({noise_outer:.3f}") fills '
          f'{100*noise_outer/image_half_width:.0f}% of the image half-width. '
          f'Consider a larger image.')
else:
    print(f'  Geometry OK: noise annulus outer edge {noise_outer:.3f}" '
          f'fills {100*noise_outer/image_half_width:.0f}% of image half-width ({image_half_width:.3f}")')

robust = config.get("robust", 0.5)
nsigma = config.get("nsigma", 3.0)

# === Step 5: Fill CASA script template ===
with open("CASA_scripts.py") as f:
    template = f.read()

filled = template.format(
    project_name    = project_name,
    output_dir      = os.path.abspath(output_dir),
    central_freq    = central_freq,
    bandwidth       = bandwidth,
    antennalist_file = antennalist_file,
    integration_time = config["integration_time"],
    pwv             = pwv,
    PA              = config["PA"],
    mask_ra         = mask_ra,
    mask_dec        = mask_dec,
    mask_semimajor  = mask_semimajor,
    mask_semiminor  = mask_semiminor,
    noise_inner     = noise_inner,
    noise_outer     = noise_outer,
    cellsize        = cellsize,
    imsize          = imsize,
    tclean_dir      = os.path.join(output_dir, "tclean"),
    base_name       = name,
    msfile          = f"{project_name}/{project_name}.alma.cycle11.{config_name.split('-')[1]}.noisy.ms",
    robust          = robust,
    nsigma          = nsigma,
)

# === Step 6: Save and run CASA ===
casa_script_path = os.path.join(output_dir, "simulate_tmp_casa.py")
with open(casa_script_path, "w") as f:
    f.write(filled)

print(f"\nDustForge single-config driver")
print(f"  Config    : {config_path}")
print(f"  Output    : {output_dir}")
print(f"  Config    : {config_name}  Band {band}  ({central_freq/1e9:.1f} GHz)")
print(f"  nsigma    : {nsigma}")
print(f"  robust    : {robust}")
print(f"\nLaunching CASA...\n")

subprocess.run([
    "casa", "--nogui", "--nologger", "--nologfile",
    "-c", casa_script_path
])
