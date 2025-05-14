import json
import os
import sys
import numpy as np
import matplotlib.pyplot as plt

# Add the 'models' directory to the Python path
sys.path.append("models")
from model_intensity_profile import gaussian_ring_profile, radial_to_image, save_image_to_fits

print('\n##### CREATING THE MODEL FITS #####\n')

# === STEP 1: Load configuration from command-line ===
if len(sys.argv) < 2:
    print("Usage: python build_model_image.py <config_file.json>")
    sys.exit(1)

config_path = sys.argv[1]
with open(config_path, "r") as f:
    config = json.load(f)

# === STEP 2: Set up 1D radial grid ===
r_max = config["disk_R90"] * 1.5
n_points = 1000
r_array = np.linspace(0, r_max, n_points)

# === STEP 3: Generate dimensionless profile ===
if config["model_type"] == "gaussian_ring":
    profile_1d = gaussian_ring_profile(
        R=r_array,
        r1=config["r1"],
        sigma1=config["sigma1"]
    )
else:
    raise ValueError("Unknown model type: {}".format(config["model_type"]))

# === STEP 4: Normalize to total flux ===
dr = r_array[1] - r_array[0]
total_raw_flux = np.sum(2 * np.pi * r_array * profile_1d * dr)
scaling_factor = config["disk_flux"] / total_raw_flux
profile_1d_scaled = profile_1d * scaling_factor

# === STEP 5: Make 2D image ===
pixel_scale = config["disk_R90"] / 128  # gives ~256 pixels across R90
image_2d = radial_to_image(
    profile_r=profile_1d_scaled,
    r_array=r_array,
    npix=512,
    pixel_scale=pixel_scale,
    inclination=config["inclination"], 
    PA=config["PA"]
)
image_2d *= config["disk_flux"] / np.sum(image_2d)

# === STEP 6: Set up output directory ===
config_name = os.path.splitext(os.path.basename(config_path))[0]
output_dir = os.path.join("output", config_name)
os.makedirs(output_dir, exist_ok=True)

# === STEP 7: Save FITS ===
fits_path = os.path.join(output_dir, "model_image.fits")
save_image_to_fits(
    image=image_2d,
    pixel_scale_au=pixel_scale,
    output_path=fits_path,
    config=config,
    npix=512
)

# === STEP 8: Save 1D profile data ===
profile_txt_path = os.path.join(output_dir, "profile_1d.txt")
np.savetxt(profile_txt_path, np.column_stack((r_array, profile_1d_scaled)), header="radius_AU intensity_Jy_per_sr")

# === STEP 9: Optional: save profile plot ===
if config.get("save_profile_plot", False):
    plot_path = os.path.join(output_dir, "profile_plot.png")
    plt.figure()
    plt.plot(r_array, profile_1d_scaled)
    plt.xlabel("Radius [AU]")
    plt.ylabel("Intensity [Jy/sr]")
    plt.title("1D Radial Profile")
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(plot_path)
    plt.close()
    print("Saved profile plot to:", plot_path)

# === STEP 10: Save a copy of the config ===
config_copy_path = os.path.join(output_dir, "config_used.json")
with open(config_copy_path, "w") as f:
    json.dump(config, f, indent=4)

print("Finished! Disk model saved to:", fits_path)
