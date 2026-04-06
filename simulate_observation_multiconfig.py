"""
simulate_observation_multiconfig.py
====================================
DustForge driver for multi-configuration ALMA synthetic observations.

Usage
-----
    python simulate_observation_multiconfig.py configs/my_proposal.json

The JSON config must include a "configurations" list, e.g.:

    "configurations": [
        {"ALMA_configuration": "C-8", "integration_time": 2448},
        {"ALMA_configuration": "C-6", "integration_time": 1804},
        {"ALMA_configuration": "C-3", "integration_time": 902}
    ]

All other top-level keys (ALMA_band, RA, Dec, disk_R90, distance, etc.)
follow the same convention as the single-config pipeline.

Workflow
--------
For each configuration:
    1. simobserve  → one noisy.ms per configuration
    2. First-pass tclean (150 iters) → RMS estimation
    3. SetWeights  → replace uniform weights with noise-based weights
Then:
    4. concat      → single merged MS
    5. Joint tclean on the concatenated MS (threshold = N × combined RMS)
    6. Export FITS + image statistics
    7. Export UV table for frank
"""

import os
import sys
import json
import subprocess
from astropy.coordinates import SkyCoord
import astropy.units as u
import numpy as np
from ALMA_specs import alma_dict

# ------------------------------------------------------------------ #
# Step 1: Read config
# ------------------------------------------------------------------ #
if len(sys.argv) < 2:
    print("Usage: python simulate_observation_multiconfig.py <config.json>")
    sys.exit(1)

config_path = sys.argv[1]
with open(config_path) as f:
    config = json.load(f)

# Validate multi-config key
if "configurations" not in config:
    print("ERROR: config must contain a 'configurations' list for the multi-config pipeline.")
    print("       Use simulate_observation.py for single-configuration runs.")
    sys.exit(1)

name       = os.path.splitext(os.path.basename(config_path))[0]
output_dir = os.path.join("output", name)
os.makedirs(output_dir, exist_ok=True)

# ------------------------------------------------------------------ #
# Step 2: Fill band-level defaults from ALMA_specs
# ------------------------------------------------------------------ #
band         = config["ALMA_band"]
specs        = alma_dict[band]
central_freq = config.get("central_freq", specs["central_freq"])
bandwidth    = config.get("bandwidth",    specs["bandwidth"])
pwv          = config.get("pwv",          specs["pwv"])   # default; can be overridden per-config

# ------------------------------------------------------------------ #
# Step 3: Build per-configuration parameter lists
# ------------------------------------------------------------------ #
configurations = config["configurations"]

project_names     = []
antennalist_files = []
integration_times = []
pwv_values        = []
ms_files          = []
cfg_labels        = []

for cfg_entry in configurations:
    cfg_name = cfg_entry["ALMA_configuration"]          # e.g. "C-8"
    int_time = cfg_entry["integration_time"]            # seconds (int or float)
    cfg_pwv  = cfg_entry.get("pwv", pwv)               # per-config PWV override

    cycle_num = cfg_name.split("-")[1]                  # "8"
    proj      = f"{config['model_type']}_B{band}_{cfg_name}_t{int_time}s"
    antlist   = f"../../alma-configuration-files/alma.cycle11.{cycle_num}.cfg"
    ms_path   = f"{proj}/{proj}.alma.cycle11.{cycle_num}.noisy.ms"

    project_names.append(proj)
    antennalist_files.append(antlist)
    integration_times.append(f"{int_time}s")
    pwv_values.append(cfg_pwv)
    ms_files.append(ms_path)
    cfg_labels.append(cfg_name)

# ------------------------------------------------------------------ #
# Step 4: Imaging geometry
# ------------------------------------------------------------------ #
# Use the highest-resolution config for the cellsize calculation
res_arcsec_values = [
    specs["resolution_arcsec"][ce["ALMA_configuration"]]
    for ce in configurations
]
min_res_arcsec = min(res_arcsec_values)  # finest beam sets cellsize

skycoord       = SkyCoord(config["RA"], config["Dec"], unit=(u.hourangle, u.deg))
mask_ra        = skycoord.ra.to_string(unit=u.hourangle, sep='hms', precision=6, pad=True)
mask_dec       = skycoord.dec.to_string(unit=u.deg, sep='dms', precision=5, alwayssign=True)
disk_size_arcsec = config["disk_R90"] / config["distance"]     # R90 in arcsec

mask_semimajor = 1.5 * disk_size_arcsec
noise_inner    = 2.5 * disk_size_arcsec
noise_outer    = 3.5 * disk_size_arcsec

# Cell size: 1/8 of the finest beam, rounded down to nearest milli-arcsec
cellsize = np.floor((min_res_arcsec / 8) * 1000) / 1000

# Image size: must contain the disk with room for the noise annulus
fov_needed = 4 * disk_size_arcsec
if fov_needed < 512 * cellsize:
    imsize = 512
elif fov_needed < 1024 * cellsize:
    imsize = 1024
else:
    imsize = 2048

robust    = config.get("robust",    0.5)
threshold = config.get("threshold", 2.0)   # factor × combined RMS

# ------------------------------------------------------------------ #
# Step 5: Fill the CASA script template
# ------------------------------------------------------------------ #
with open("CASA_scripts_multiconfig.py") as f:
    template = f.read()

# Helper: render Python list literals for CASA (which runs Python 3)
def py_list(lst):
    return repr(lst)

filled = template.format(
    output_dir        = os.path.abspath(output_dir),
    project_names     = py_list(project_names),
    antennalist_files = py_list(antennalist_files),
    integration_times = py_list(integration_times),
    pwv_values        = py_list(pwv_values),
    central_freq      = central_freq,
    bandwidth         = bandwidth,
    ms_files          = py_list(ms_files),
    cfg_labels        = py_list(cfg_labels),
    PA                = config["PA"],
    mask_ra           = mask_ra,
    mask_dec          = mask_dec,
    mask_semimajor    = mask_semimajor,
    mask_semiminor    = mask_semimajor * np.cos(np.radians(config["inclination"])),
    noise_inner       = noise_inner,
    noise_outer       = noise_outer,
    imsize            = imsize,
    cellsize          = cellsize,
    base_name         = name,
    robust            = robust,
    threshold         = threshold,
)

# ------------------------------------------------------------------ #
# Step 6: Write and execute the CASA script
# ------------------------------------------------------------------ #
casa_script_path = os.path.join(output_dir, "simulate_multiconfig_casa.py")
with open(casa_script_path, "w") as f:
    f.write(filled)

print(f"\nDustForge multi-config driver")
print(f"  Config : {config_path}")
print(f"  Output : {output_dir}")
print(f"  Configs: {cfg_labels}")
print(f"  Band   : {band}  ({central_freq/1e9:.1f} GHz)")
print(f"  Cellsize: {cellsize:.3f} arcsec  imsize: {imsize}")
print(f"  Mask semimajor: {mask_semimajor:.2f} arcsec")
print(f"\nLaunching CASA...\n")

subprocess.run([
    "casa", "--nogui", "--nologger", "--nologfile",
    "-c", casa_script_path
])
