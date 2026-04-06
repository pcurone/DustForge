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

Workflow
--------
For each configuration:
    1. simobserve  -> one noisy MS
    2. Full tclean with nsigma stopping -> per-config image
    3. Measure RMS from that image -> SetWeights on the MS
Then:
    4. concat all (now correctly weighted) MSes
    5. Joint tclean with nsigma on the concatenated MS
    6. Export FITS + stats + UV table for frank
"""

import os
import sys
import json
import re
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

if "configurations" not in config:
    print("ERROR: config must contain a 'configurations' list.")
    print("       Use simulate_observation.py for single-configuration runs.")
    sys.exit(1)

name       = os.path.splitext(os.path.basename(config_path))[0]
output_dir = os.path.join("output", name)
os.makedirs(output_dir, exist_ok=True)

# ------------------------------------------------------------------ #
# Step 2: Band-level defaults
# ------------------------------------------------------------------ #
band         = config["ALMA_band"]
specs        = alma_dict[band]
central_freq = config.get("central_freq", specs["central_freq"])
bandwidth    = config.get("bandwidth",    specs["bandwidth"])
pwv_default  = config.get("pwv",          specs["pwv"])

# ------------------------------------------------------------------ #
# Step 3: Per-configuration parameter lists
# ------------------------------------------------------------------ #
configurations    = config["configurations"]
project_names     = []
antennalist_files = []
integration_times = []
pwv_values        = []
ms_files          = []
cfg_labels        = []

for cfg_entry in configurations:
    cfg_name  = cfg_entry["ALMA_configuration"]
    int_time  = cfg_entry["integration_time"]
    cfg_pwv   = cfg_entry.get("pwv", pwv_default)
    cycle_num = cfg_name.split("-")[1]
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
res_arcsec_values = [
    specs["resolution_arcsec"][ce["ALMA_configuration"]]
    for ce in configurations
]
min_res_arcsec = min(res_arcsec_values)
# The most extended config has the finest (smallest) beam
extended_cfg_idx  = res_arcsec_values.index(min_res_arcsec)
extended_cfg_name = cfg_labels[extended_cfg_idx]

skycoord         = SkyCoord(config["RA"], config["Dec"], unit=(u.hourangle, u.deg))
mask_ra          = skycoord.ra.to_string(unit=u.hourangle, sep='hms', precision=6, pad=True)
mask_dec         = skycoord.dec.to_string(unit=u.deg, sep='dms', precision=5, alwayssign=True)
disk_size_arcsec = config["disk_R90"] / config["distance"]

mask_semimajor = 1.5 * disk_size_arcsec
mask_semiminor = mask_semimajor * np.cos(np.radians(config["inclination"]))
noise_inner    = 2.5 * disk_size_arcsec
noise_outer    = 3.5 * disk_size_arcsec

print(f"  Extended config (finest beam): {extended_cfg_name}  "
      f"({min_res_arcsec}\" beam) — cellsize/imsize apply here and on joint tclean")

# Check for inactive underscore-prefixed overrides
for k in ["_cellsize", "_imsize"]:
    if k in config:
        print(f"  NOTE: found \"{k}\" in config — INACTIVE. "
              f"Rename to \"{k[1:]}\" to activate.")

# cellsize and imsize apply to the extended config and the joint tclean only.
# Compact configs always use their own auto-computed values.
if "cellsize" in config:
    cellsize = config["cellsize"]
    print(f"  cellsize (extended + joint): MANUAL OVERRIDE = {cellsize:.4f} arcsec")
else:
    cellsize = max(np.floor((min_res_arcsec / 8) * 1000) / 1000, 0.002)
    print(f"  cellsize (extended + joint): AUTO = {cellsize:.4f} arcsec  (beam/8 = {min_res_arcsec/8:.4f})")

if "imsize" in config:
    imsize = config["imsize"]
    print(f"  imsize (extended + joint): MANUAL OVERRIDE = {imsize}")
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
    print(f"  imsize (extended + joint): AUTO = {imsize}  "
          f"(FOV needed: {fov_needed:.3f}\"  covers: {imsize*cellsize:.3f}\")")


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

# ------------------------------------------------------------------ #
# Step 5: Fill CASA template using <<<token>>> substitution
# ------------------------------------------------------------------ #
with open("CASA_scripts_multiconfig.py") as f:
    script = f.read()

substitutions = {
    "<<<output_dir>>>":        os.path.abspath(output_dir),
    "<<<project_names>>>":     repr(project_names),
    "<<<antennalist_files>>>": repr(antennalist_files),
    "<<<integration_times>>>": repr(integration_times),
    "<<<pwv_values>>>":        repr(pwv_values),
    "<<<central_freq>>>":      str(int(central_freq)),
    "<<<bandwidth>>>":         str(int(bandwidth)),
    "<<<ms_files>>>":          repr(ms_files),
    "<<<cfg_labels>>>":        repr(cfg_labels),
    "<<<PA>>>":                f"{config['PA']:.1f}",
    "<<<mask_ra>>>":           mask_ra,
    "<<<mask_dec>>>":          mask_dec,
    "<<<mask_semimajor>>>":    f"{mask_semimajor:.3f}arcsec",
    "<<<mask_semiminor>>>":    f"{mask_semiminor:.3f}arcsec",
    "<<<noise_inner>>>":       f"{noise_inner:.3f}",
    "<<<noise_outer>>>":       f"{noise_outer:.3f}",
    "<<<res_arcsec_per_cfg>>>": repr(res_arcsec_values),
    "<<<extended_cfg>>>":      repr(extended_cfg_name),
    "<<<imsize>>>":            str(imsize),
    "<<<cellsize>>>":          f"{cellsize:.4f}",
    "<<<base_name>>>":         name,
    "<<<robust>>>":            str(robust),
    "<<<nsigma>>>":            str(nsigma),
}

for token, value in substitutions.items():
    script = script.replace(token, value)

# Safety check
remaining = re.findall(r'<<<\w+>>>', script)
if remaining:
    print(f"WARNING: unreplaced template tokens: {remaining}")

# ------------------------------------------------------------------ #
# Step 6: Write and run CASA script
# ------------------------------------------------------------------ #
casa_script_path = os.path.join(output_dir, "simulate_multiconfig_casa.py")
with open(casa_script_path, "w") as f:
    f.write(script)

print(f"\nDustForge multi-config driver")
print(f"  Config     : {config_path}")
print(f"  Output     : {output_dir}")
print(f"  Configs    : {cfg_labels}")
print(f"  Band       : {band}  ({central_freq/1e9:.1f} GHz)")
print(f"  Cellsize (extended+joint): {cellsize:.4f} arcsec  imsize: {imsize}")
print(f"  Compact configs          : auto cellsize/imsize per config")
print(f"  nsigma     : {nsigma}")
print(f"  robust     : {robust}")
print(f"  Mask smaj  : {mask_semimajor:.3f} arcsec  smin: {mask_semiminor:.3f} arcsec")
print(f"  Noise annulus: {noise_inner:.3f} – {noise_outer:.3f} arcsec")
print(f"\nLaunching CASA...\n")

subprocess.run([
    "casa", "--nogui", "--nologger", "--nologfile",
    "-c", casa_script_path
])
