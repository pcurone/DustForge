import os
import sys
import json
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from gofish import imagecube

run_frank = True

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
incl = config["inclination"]
PA = config["PA"]

# === Step 2: Compute CLEAN intensity radial profile ===

# load data
data_fits = os.path.join(output_dir, "tclean", f"{name}_robust0.5_clean.image.fits") 
dhdu = fits.open(data_fits)
dimg, hd = np.squeeze(dhdu[0].data), dhdu[0].header
bmaj, bmin, bPA = 3600 * hd['BMAJ'], 3600 * hd['BMIN'], hd['BPA']
beam_area_sr = (np.pi * bmaj * bmin / (4 * np.log(2))) / (3600 * 180 / np.pi)**2

# Obtain the CLEAN profile using the imagecube function from gofish
FOV_gofish = 5     # Clip the image cube down to a specific field-of-view spanning a range ``FOV``, where ``FOV`` is in [arcsec]
cube = imagecube(data_fits, FOV=FOV_gofish)
r_clean, I_clean, dI_clean = cube.radial_profile(x0=0.0, y0=0.0, inc=incl, PA=PA)
r_clean_au = r_clean * distance_pc
I_clean_jy_sr = I_clean / beam_area_sr
dI_clean_jy_sr = dI_clean / beam_area_sr

# Save CLEAN radial profile
profile_path = os.path.join(output_dir, "CLEAN_profile.txt")
np.savetxt(profile_path, np.column_stack([r_clean_au, I_clean_jy_sr, dI_clean_jy_sr]),
            header=" Radius(au)\tIntensity(Jy/sr)\td_Intensity(Jy/sr)")


# === Step 3: Compute recovery score for the CLEAN intensity radial profile ===

def compute_recovery_score(model_r, model_I, data_r, data_I, data_err, save_profiles=True):
    """
    Compute two recovery scores comparing the model and recovered (CLEAN) radial intensity profiles,
    taking into account the uncertainty on the CLEAN profile via inverse-variance weighting.

    - Weighted L1 score (normalized weighted mean absolute error):
      Emphasizes overall similarity; robust to occasional mismatches and local deviations.

    - Weighted L2 score (normalized weighted root mean square error):
      Penalizes large errors more strongly; sensitive to sharp differences like missed peaks or gaps.

    Both scores range from:
    - 1: perfect recovery
    - 0: no recovery
    - <0: very poor recovery (model deviates strongly from CLEAN profile)
    """

    # Interpolate CLEAN profile to model radii
    data_I_interp = np.interp(model_r, data_r, data_I, left=np.nan, right=np.nan)
    data_err_interp = np.interp(model_r, data_r, data_err, left=np.nan, right=np.nan)

    # Mask NaNs
    mask = ~np.isnan(data_I_interp) & ~np.isnan(data_err_interp)
    model_I = model_I[mask]
    data_I_interp = data_I_interp[mask]
    data_err_interp = data_err_interp[mask]

    # Weighted L2 score
    l2_numerator = np.sum(((model_I - data_I_interp)**2) / data_err_interp**2)
    l2_denominator = np.sum((model_I**2) / data_err_interp**2)
    l2_score = 1 - np.sqrt(l2_numerator / l2_denominator)

    # Weighted L1 score
    l1_numerator = np.sum(np.abs(model_I - data_I_interp) / data_err_interp)
    l1_denominator = np.sum(np.abs(model_I) / data_err_interp)
    l1_score = 1 - (l1_numerator / l1_denominator)

    # Save interpolated comparison arrays for inspection
    if save_profiles:
        comparison_data = np.column_stack([model_r[mask], model_I, data_I_interp, data_err_interp])
        comparison_path = os.path.join(output_dir, "interpolated_profiles_model_cleanI.txt")
        np.savetxt(comparison_path, comparison_data,
                  header="# Radius [AU]    Model_I [Jy/sr]    data_I_interp [Jy/sr]    data_err_interp [Jy/sr]",
                  fmt="%.6e")

    return round(l1_score, 4), round(l2_score, 4)
    
# Load model 1D profile
model_profile_path = os.path.join(output_dir, "profile_model.txt")
model_data = np.loadtxt(model_profile_path)
r_model_au, I_model_jy_sr = model_data[:, 0], model_data[:, 1]

# Compute recovery scores
l1_score, l2_score = compute_recovery_score(r_model_au, I_model_jy_sr,
                                            r_clean_au, I_clean_jy_sr, dI_clean_jy_sr)


# === Step 4: Plot the model and CLEAN profiles ===

plot_path = os.path.join(output_dir, "profile_comparison_model_CLEAN.png")
plt.figure()

plt.plot(r_model_au, I_model_jy_sr, label="Input Model", lw=2)
plt.plot(r_clean_au, I_clean_jy_sr, label="CLEAN Profile", lw=2)
plt.fill_between(r_clean_au,
                 I_clean_jy_sr - dI_clean_jy_sr,
                 I_clean_jy_sr + dI_clean_jy_sr, 
                 color="tab:orange", alpha=0.3, label="CLEAN ±1σ")

# Add labels and grid
plt.xlim(0, 300)
plt.xlabel("Radius [AU]")
plt.ylabel("Intensity [Jy/sr]")
plt.title("Model vs CLEAN Radial Profile")
plt.grid(True)
plt.legend()

# Add recovery scores as text inside the plot
text_str = f"L1 recovery: {l1_score:.3f}\nL2 recovery: {l2_score:.3f}"
plt.text(0.95, 0.20, text_str,
         transform=plt.gca().transAxes,
         ha='right',
         fontsize=10,
         verticalalignment='top',
         bbox=dict(facecolor='white', edgecolor='gray', alpha=0.8))

# Finalize and save
plt.tight_layout()
plt.savefig(plot_path)
plt.close()
print("Saved profile comparison plot to:", plot_path)

# === Step 5: Save recovery scores ===
score_path = os.path.join(output_dir, "recovery_scores_CLEAN.txt")
with open(score_path, "w") as f:
    f.write(f"L1 recovery score (CLEAN): {l1_score}\n")
    f.write(f"L2 recovery score (CLEAN): {l2_score}\n")



# === Step 6: Optional FRANK module ===
if run_frank:

  import frank
  from frank.radial_fitters import FrankFitter
  from frank.geometry import FixedGeometry
  from frank.io import save_fit, load_sol
  from frank.make_figs import make_full_fig
  frank.enable_logging()


  print(' ')
  print('##### RUNNING FRANK #####')
  print(' ')

  # load the visibility data
  frank_dir = os.path.join(output_dir, "frank_fit")
  dat = np.load(os.path.join(frank_dir, "uvtable_frank.npz"))
  u, v, vis, wgt = dat['u'], dat['v'], dat['Vis'], dat['Wgt']

  # set the disk viewing geometry
  geom = FixedGeometry(incl, PA, dRA=0.0, dDec=0.0)

  # configure the fitting code setup
  FF = FrankFitter(Rmax=config["disk_R90"]/config["distance"]*2, geometry=geom, N=300, alpha=1.30, weights_smooth=0.01, method='LogNormal')
  
  # fit the visibilities
  sol = FF.fit(u, v, vis, wgt)
  make_full_fig(u, v, vis, wgt, sol, bin_widths=[1e4, 5e4], save_prefix=f'{frank_dir}/summary')

  # read results
  r_frank_arcsec = sol.r
  r_frank_au = r_frank_arcsec * distance_pc  # from arcsec to AU
  I_frank = sol.I

  # save the fit
  np.savetxt(os.path.join(frank_dir, f"frank_profile.txt"),
              np.column_stack([r_frank_arcsec, r_frank_au, I_frank]),
              header="# R(arcsec)\tR(au)\tIntensity(Jy/sr)", fmt="%.6e")

  l1_frank, l2_frank = compute_recovery_score(r_model_au, I_model_jy_sr,
                                              r_frank_au, I_frank, np.ones_like(I_frank))

  frank_plot_path = os.path.join(frank_dir, "profile_comparison_model_CLEAN_frank.png")
  plt.figure()
  plt.plot(r_model_au, I_model_jy_sr, label="Input Model", lw=2)
  plt.plot(r_clean_au, I_clean_jy_sr, label="CLEAN Profile", lw=2)
  plt.fill_between(r_clean_au, I_clean_jy_sr - dI_clean_jy_sr, I_clean_jy_sr + dI_clean_jy_sr,
                    color="tab:orange", alpha=0.3)
  plt.plot(r_frank_au, I_frank, label="frank Profile", lw=2, ls="--")
  plt.xlabel("Radius [AU]")
  plt.ylabel("Intensity [Jy/sr]")
  plt.title("Model vs CLEAN vs frank Radial Profile")
  plt.grid(True)
  plt.legend()
  plt.text(0.95, 0.25,
            f"L1 CLEAN: {l1_score:.3f}\nL2 CLEAN: {l2_score:.3f}\nL1 frank: {l1_frank:.3f}\nL2 frank: {l2_frank:.3f}",
            transform=plt.gca().transAxes,
            ha='right', fontsize=10, verticalalignment='top',
            bbox=dict(facecolor='white', edgecolor='gray', alpha=0.8))
  plt.tight_layout()
  plt.savefig(frank_plot_path)
  plt.close()
  print("Saved model vs CLEAN vs frank comparison to:", frank_plot_path)

  with open(os.path.join(frank_dir, "recovery_scores_frank.txt"), "w") as f:
      f.write("# Columns: L1_recovery_frank\tL2_recovery_frank\n")
      f.write(f"{l1_frank:.4f}\t{l2_frank:.4f}\n")
