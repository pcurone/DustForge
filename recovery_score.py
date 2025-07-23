import os
import sys
import json
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from gofish import imagecube

# === Step 1: Read config ===
if len(sys.argv) < 2:
    print("Usage: python simulate_observation.py <config.json>")
    sys.exit(1)

config_path = sys.argv[1]
with open(config_path) as f:
    config = json.load(f)

name = os.path.splitext(os.path.basename(config_path))[0]
output_dir = os.path.join("output", name)
distance_pc = config["distance"]
incl = config["incl"]
PA = config["PA"]

# === Step 2: Compute CLEAN intensity radial profile ===

# load data
data_fits = os.path.join(output_dir, "tclean", f"{name}_robust0.5.clean.image.fits") 
dhdu = fits.open(data_fits)
dimg, hd = np.squeeze(dhdu[0].data), dhdu[0].header
bmaj, bmin, bPA = 3600 * hd['BMAJ'], 3600 * hd['BMIN'], hd['BPA']
beam_area_sr = (np.pi * bmaj * bmin / (4 * np.log(2))) / (3600 * 180 / np.pi)**2

# Obtain the CLEAN profile using the imagecube function from gofish
cube = imagecube(data_fits, FOV=FOV_gofish)
r_clean, I_clean, dI_clean = cube.radial_profile(x0=0.0, y0=0.0, inc=inclination, PA=PA)
r_clean_au = r_clean * distance_pc
I_clean_jy_sr = I_clean / beam_area_sr
dI_clean_jy_sr = dI_clean / beam_area_sr

# Save CLEAN radial profile
profile_path = os.path.join(output_dir, "CLEAN_profile.txt")
np.savetxt(profile_path, np.column_stack([r_clean_au, I_clean_jy_sr, dI_clean_jy_sr]),
            header="# Radius(au)\tIntensity(Jy/sr)\td_Intensity(Jy/sr)")


# === Step 3: Compute recovery score for the CLEAN intensity radial profile ===

def compute_recovery_score(model_r, model_I, clean_r, clean_I):
    """
    Compute two normalized recovery scores comparing the model and recovered radial intensity profiles.

    - L1 score (based on normalized mean absolute error):
      Emphasizes overall similarity; robust to occasional mismatches and large deviations.

    - L2 score (based on normalized root mean square error):
      Penalizes large errors more strongly; sensitive to sharp differences like missed peaks or gaps.

    Both scores range from 0 (no recovery) to 1 (perfect recovery).
    """
    # Interpolate CLEAN profile to model radii
    clean_I_interp = np.interp(model_r, clean_r, clean_I, left=np.nan, right=np.nan)

    # L2 recovery score
    diff_l2 = model_I - clean_I
    l2_score = 1 - np.sqrt(np.mean(diff_l2**2))

    # L1 recovery score
    diff_l1 = np.abs(model_norm - clean_norm)
    l1_score = 1 - np.mean(diff_l1)

    return round(l1_score, 4), round(l2_score, 4
    
# Load model 1D profile
model_profile_path = os.path.join(output_dir, "profile_1D.txt")
model_data = np.loadtxt(model_profile_path)
r_model_au, I_model_jy_sr = model_data[:, 0], model_data[:, 1]

# Compute recovery scores
l1_score, l2_score = compute_recovery_score(r_model_au, I_model_jy_sr,
                                            r_clean_au, I_clean_jy_sr)


plot_path = os.path.join(output_dir, "profile_comparison.png")
plt.figure()
plt.plot(r_model_au, I_model_jysr, label="Input Model", lw=2)
plt.plot(r_clean_au, I_clean_jysr, label="CLEAN Profile", lw=2)
plt.fill_between(r_clean_au,
                    I_clean_jysr - dI_clean_jysr,
                    I_clean_jysr + dI_clean_jysr,
                    color="gray", alpha=0.3, label="CLEAN ±1σ")

plt.xlabel("Radius [AU]")
plt.ylabel("Intensity [Jy/sr]")
plt.title("Model vs CLEAN Radial Profile")
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.savefig(plot_path)
plt.close()
print("Saved profile comparison plot to:", plot_path)

# Save recovery scores
score_path = os.path.join(output_dir, "recovery_scores.txt")
with open(score_path, "w") as f:
    f.write(f"L1 recovery score (CLEAN): {l1_score}\n")
    f.write(f"L2 recovery score (CLEAN): {l2_score}\n")






















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

cellsize = np.floor((res_arcsec / 6) * 1000) / 1000  # round down to nearest milliarcsec

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
