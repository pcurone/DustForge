"""
recovery_score.py
=================
Computes radial profile recovery scores comparing the input model
against the CLEAN image (and optionally frank) for either:
  - Single-config runs  (reads tclean/{name}_robust{r}_clean.image.fits)
  - Multi-config runs   (reads tclean/{name}_multiconfig_robust{r}_clean.image.fits)

The script auto-detects which image exists.

Usage
-----
    python recovery_score.py configs/my_config.json
"""

import os
import sys
import json
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from gofish import imagecube

run_frank = True

# ------------------------------------------------------------------ #
# Step 1: Read config
# ------------------------------------------------------------------ #
if len(sys.argv) < 2:
    print("Usage: python recovery_score.py <config.json>")
    sys.exit(1)

config_path = sys.argv[1]
with open(config_path) as f:
    config = json.load(f)

name        = os.path.splitext(os.path.basename(config_path))[0]
output_dir  = os.path.join("output", name)
distance_pc = config["distance"]
incl        = config["inclination"]
PA          = config["PA"]
robust      = config.get("robust", 0.5)

# ------------------------------------------------------------------ #
# Step 2: Locate the CLEAN image (single- or multi-config)
# ------------------------------------------------------------------ #
tclean_dir = os.path.join(output_dir, "tclean")

single_fits = os.path.join(tclean_dir, f"{name}_robust{robust}_clean.image.fits")
multi_fits  = os.path.join(tclean_dir, f"{name}_multiconfig_robust{robust}_clean.image.fits")

if os.path.exists(multi_fits):
    data_fits = multi_fits
    mode_label = "multiconfig"
    print(f"[recovery_score] Multi-config image found: {data_fits}")
elif os.path.exists(single_fits):
    data_fits = single_fits
    mode_label = "single"
    print(f"[recovery_score] Single-config image found: {data_fits}")
else:
    print(f"ERROR: No CLEAN FITS image found in {tclean_dir}")
    print(f"  Looked for:\n    {single_fits}\n    {multi_fits}")
    sys.exit(1)

# ------------------------------------------------------------------ #
# Step 3: Compute CLEAN radial profile
# ------------------------------------------------------------------ #
dhdu = fits.open(data_fits)
dimg, hd = np.squeeze(dhdu[0].data), dhdu[0].header
bmaj, bmin, bPA = 3600 * hd['BMAJ'], 3600 * hd['BMIN'], hd['BPA']
beam_area_sr = (np.pi * bmaj * bmin / (4 * np.log(2))) / (3600 * 180 / np.pi)**2

FOV_gofish = 5   # arcsec — clip to disk region
cube = imagecube(data_fits, FOV=FOV_gofish)
r_clean, I_clean, dI_clean = cube.radial_profile(x0=0.0, y0=0.0, inc=incl, PA=PA)
r_clean_au    = r_clean * distance_pc
I_clean_jy_sr = I_clean  / beam_area_sr
dI_clean_jy_sr = dI_clean / beam_area_sr

profile_path = os.path.join(output_dir, f"CLEAN_profile_{mode_label}.txt")
np.savetxt(profile_path,
           np.column_stack([r_clean_au, I_clean_jy_sr, dI_clean_jy_sr]),
           header="Radius(au)\tIntensity(Jy/sr)\td_Intensity(Jy/sr)")

# ------------------------------------------------------------------ #
# Step 4: Recovery score
# ------------------------------------------------------------------ #
def compute_recovery_score(model_r, model_I, data_r, data_I, data_err,
                            output_dir=output_dir, tag=""):
    """
    Weighted L1 and L2 recovery scores.

    L1 score = 1 - sum(|model - data| / err) / sum(|model| / err)
    L2 score = 1 - sqrt(sum((model-data)^2/err^2) / sum(model^2/err^2))

    Both equal 1 for perfect recovery, 0 for no correlation.
    Scores < 0 indicate the model is a worse predictor than zero.
    """
    data_I_interp   = np.interp(model_r, data_r, data_I,   left=np.nan, right=np.nan)
    data_err_interp = np.interp(model_r, data_r, data_err, left=np.nan, right=np.nan)

    mask = ~np.isnan(data_I_interp) & ~np.isnan(data_err_interp) & (data_err_interp > 0)
    model_I_m        = model_I[mask]
    data_I_interp_m  = data_I_interp[mask]
    data_err_interp_m = data_err_interp[mask]

    l2_num   = np.sum(((model_I_m - data_I_interp_m)**2) / data_err_interp_m**2)
    l2_denom = np.sum((model_I_m**2) / data_err_interp_m**2)
    l2_score = 1 - np.sqrt(l2_num / l2_denom)

    l1_num   = np.sum(np.abs(model_I_m - data_I_interp_m) / data_err_interp_m)
    l1_denom = np.sum(np.abs(model_I_m) / data_err_interp_m)
    l1_score = 1 - (l1_num / l1_denom)

    # Save interpolated comparison
    comparison_path = os.path.join(output_dir, f"interpolated_profiles_{tag}.txt")
    np.savetxt(comparison_path,
               np.column_stack([model_r[mask], model_I_m, data_I_interp_m, data_err_interp_m]),
               header="Radius[AU]  Model_I[Jy/sr]  data_I[Jy/sr]  data_err[Jy/sr]",
               fmt="%.6e")

    return round(l1_score, 4), round(l2_score, 4)


model_data    = np.loadtxt(os.path.join(output_dir, "profile_model.txt"))
r_model_au    = model_data[:, 0]
I_model_jy_sr = model_data[:, 1]

l1_score, l2_score = compute_recovery_score(
    r_model_au, I_model_jy_sr,
    r_clean_au, I_clean_jy_sr, dI_clean_jy_sr,
    tag=f"model_cleanI_{mode_label}"
)

# ------------------------------------------------------------------ #
# Step 5: Plot model vs CLEAN
# ------------------------------------------------------------------ #
plot_path = os.path.join(output_dir, f"profile_comparison_model_CLEAN_{mode_label}.png")
fig, ax = plt.subplots()
ax.plot(r_model_au, I_model_jy_sr, label="Input Model", lw=2)
ax.plot(r_clean_au, I_clean_jy_sr, label="CLEAN Profile", lw=2)
ax.fill_between(r_clean_au,
                I_clean_jy_sr - dI_clean_jy_sr,
                I_clean_jy_sr + dI_clean_jy_sr,
                color="tab:orange", alpha=0.3, label="CLEAN ±1σ")
ax.set_xlim(0, 300)
ax.set_xlabel("Radius [AU]")
ax.set_ylabel("Intensity [Jy/sr]")
ax.set_title(f"Model vs CLEAN — {mode_label}")
ax.grid(True)
ax.legend()
ax.text(0.95, 0.20, f"L1 recovery: {l1_score:.3f}\nL2 recovery: {l2_score:.3f}",
        transform=ax.transAxes, ha='right', fontsize=10, verticalalignment='top',
        bbox=dict(facecolor='white', edgecolor='gray', alpha=0.8))
plt.tight_layout()
plt.savefig(plot_path, dpi=150)
plt.close()
print(f"Saved profile comparison: {plot_path}")

# ------------------------------------------------------------------ #
# Step 6: Save CLEAN recovery scores
# ------------------------------------------------------------------ #
score_path = os.path.join(output_dir, f"recovery_scores_CLEAN_{mode_label}.txt")
with open(score_path, "w") as f:
    f.write(f"Mode: {mode_label}\n")
    f.write(f"L1 recovery score (CLEAN): {l1_score}\n")
    f.write(f"L2 recovery score (CLEAN): {l2_score}\n")

# ------------------------------------------------------------------ #
# Step 7: Optional frank module
# ------------------------------------------------------------------ #
if run_frank:
    import frank
    from frank.radial_fitters import FrankFitter
    from frank.geometry import FixedGeometry
    from frank.make_figs import make_full_fig
    frank.enable_logging()

    print('\n##### RUNNING FRANK #####\n')

    frank_dir = os.path.join(output_dir, "frank_fit")
    os.makedirs(frank_dir, exist_ok=True)

    # Prefer the multi-config UV table if present
    uvtable_multi  = os.path.join(frank_dir, "uvtable_frank_multiconfig.npz")
    uvtable_single = os.path.join(frank_dir, "uvtable_frank.npz")

    if os.path.exists(uvtable_multi):
        uvtable_path = uvtable_multi
        print(f"  Using multi-config UV table: {uvtable_path}")
    elif os.path.exists(uvtable_single):
        uvtable_path = uvtable_single
        print(f"  Using single-config UV table: {uvtable_path}")
    else:
        print("  WARNING: No UV table found — skipping frank.")
        run_frank = False

if run_frank:
    dat = np.load(uvtable_path)
    u_frank, v_frank = dat['u'], dat['v']
    vis_frank, wgt_frank = dat['Vis'], dat['Wgt']

    geom = FixedGeometry(incl, PA, dRA=0.0, dDec=0.0)
    FF   = FrankFitter(
        Rmax=config["disk_R90"] / config["distance"] * 2,
        geometry=geom,
        N=300,
        alpha=1.30,
        weights_smooth=0.01,
        method='LogNormal'
    )
    sol = FF.fit(u_frank, v_frank, vis_frank, wgt_frank)
    make_full_fig(u_frank, v_frank, vis_frank, wgt_frank, sol,
                  bin_widths=[1e4, 5e4],
                  save_prefix=os.path.join(frank_dir, "summary"))

    r_frank_arcsec = sol.r
    r_frank_au     = r_frank_arcsec * distance_pc
    I_frank        = sol.I

    np.savetxt(os.path.join(frank_dir, "frank_profile.txt"),
               np.column_stack([r_frank_arcsec, r_frank_au, I_frank]),
               header="R(arcsec)\tR(au)\tIntensity(Jy/sr)", fmt="%.6e")

    l1_frank, l2_frank = compute_recovery_score(
        r_model_au, I_model_jy_sr,
        r_frank_au, I_frank, np.ones_like(I_frank),
        tag="model_frank"
    )

    frank_plot_path = os.path.join(frank_dir, f"profile_comparison_model_CLEAN_frank_{mode_label}.png")
    fig, ax = plt.subplots()
    ax.plot(r_model_au, I_model_jy_sr, label="Input Model", lw=2)
    ax.plot(r_clean_au, I_clean_jy_sr, label="CLEAN Profile", lw=2)
    ax.fill_between(r_clean_au, I_clean_jy_sr - dI_clean_jy_sr, I_clean_jy_sr + dI_clean_jy_sr,
                    color="tab:orange", alpha=0.3)
    ax.plot(r_frank_au, I_frank, label="frank Profile", lw=2, ls="--")
    ax.set_xlabel("Radius [AU]")
    ax.set_ylabel("Intensity [Jy/sr]")
    ax.set_title(f"Model vs CLEAN vs frank — {mode_label}")
    ax.grid(True)
    ax.legend()
    ax.text(0.95, 0.28,
            f"L1 CLEAN: {l1_score:.3f}\nL2 CLEAN: {l2_score:.3f}\n"
            f"L1 frank: {l1_frank:.3f}\nL2 frank: {l2_frank:.3f}",
            transform=ax.transAxes, ha='right', fontsize=10, verticalalignment='top',
            bbox=dict(facecolor='white', edgecolor='gray', alpha=0.8))
    plt.tight_layout()
    plt.savefig(frank_plot_path, dpi=150)
    plt.close()
    print(f"Saved frank comparison: {frank_plot_path}")

    with open(os.path.join(frank_dir, "recovery_scores_frank.txt"), "w") as f:
        f.write(f"Mode: {mode_label}\n")
        f.write(f"L1 recovery (frank): {l1_frank:.4f}\n")
        f.write(f"L2 recovery (frank): {l2_frank:.4f}\n")
