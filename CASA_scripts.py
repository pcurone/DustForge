import numpy as np
import os

print('\n##### RUNNING SIMOBSERVE #####\n')

os.chdir("{output_dir}")

simobserve(
    project="{project_name}",
    skymodel="model_image.fits",
    indirection="",
    incell="",
    inbright="",
    incenter="{central_freq}Hz",
    inwidth="{bandwidth}Hz",
    antennalist="../../alma-configuration-files/{antennalist_file}",
    totaltime="{integration_time}s",
    user_pwv={pwv},
    overwrite=True
)

print('\n##### RUNNING TCLEAN #####\n')

# Mask and noise annulus geometry
mask = f"ellipse[[{mask_ra},{mask_dec}], [{mask_semimajor:.3f}arcsec, {mask_semiminor:.3f}arcsec], {PA:.1f}deg]"
noise_annulus = f"annulus[[{mask_ra}, {mask_dec}],['{noise_inner:.3f}arcsec', '{noise_outer:.3f}arcsec']]"

imsize  = {imsize}
nsigma  = {nsigma}
robust  = {robust}

os.makedirs('tclean', exist_ok=True)
image_base  = "tclean/{base_name}"
image_final = f"{{image_base}}_robust{robust}_clean"

for ext in [".image", ".mask", ".model", ".pb", ".psf", ".residual", ".sumwt", "weight"]:
    os.system(f"rm -rf {{image_final}}{{ext}}")

# ----------------------------------------------------------
# Single tclean call.
# threshold='0.0Jy' disables the fixed threshold entirely.
# nsigma lets tclean compute an N-sigma stopping criterion
# at each major cycle from the residual image noise.
# fastnoise=False uses MAD-based noise estimation, which is
# robust against bright sources biasing the RMS upward.
# ----------------------------------------------------------
print(f"tclean: nsigma={{nsigma}}  robust={robust}  fastnoise=False")

tclean(
    vis="{msfile}",
    imagename=image_final,
    datacolumn='data',
    specmode='mfs',
    deconvolver='multiscale',
    scales=[0, 8, 15, 30],
    gridder='standard',
    weighting='briggs',
    robust=robust,
    gain=0.1,
    threshold='0.0Jy',
    nsigma=nsigma,
    fastnoise=False,
    cycleniter=300,
    cyclefactor=1,
    imsize=imsize,
    cell=f"{cellsize:.3f}arcsec",
    niter=50000,
    mask=mask,
    interactive=False,
    savemodel='none',
    nterms=1
)

# Export FITS
exportfits(
    imagename=image_final + '.image',
    fitsimage=image_final + '.image.fits',
    overwrite=True,
    history=True
)

# ----------------------------------------------------------
# Measure RMS from the final image and run SetWeights.
# Even for a single MS this is useful: it leaves the MS with
# correct noise-based weights so it can be concatenated with
# other configurations later without re-processing.
# ----------------------------------------------------------
rms_jy  = imstat(imagename=image_final + '.image', region=noise_annulus)['rms'][0]
rms_mJy = rms_jy * 1e3
print(f"Measured RMS: {{rms_mJy:.3f}} mJy  --> running SetWeights")

def SetWeights(MS_filename, rms_jy):
    tb.open(MS_filename, nomodify=False)
    weight = tb.getcol("WEIGHT")
    sigma  = tb.getcol("SIGMA")
    if weight.ndim == 2:
        nvis = weight.shape[1]
    else:
        nvis = weight.shape[0]
    weight[:] = 1.0 / (2.0 * nvis * rms_jy**2)
    sigma      = 1.0 / np.sqrt(weight)
    tb.putcol("WEIGHT", weight)
    tb.putcol("SIGMA",  sigma)
    tb.flush()
    tb.close()
    print(f"SetWeights: {{MS_filename}}  rms={{rms_jy*1e3:.3f}} mJy  weight={{weight.flat[0]:.4e}}")

SetWeights("{msfile}", rms_jy)

# Beam and source statistics
headerlist     = imhead(image_final + '.image', mode='list')
beammajor      = headerlist['beammajor']['value']
beamminor      = headerlist['beamminor']['value']
beampa         = headerlist['beampa']['value']
target_stats   = imstat(imagename=image_final + '.image', region=mask)
target_flux    = target_stats['flux'][0]
peak_intensity = target_stats['max'][0]
SNR            = peak_intensity / rms_jy

with open(f'Info_image_{base_name}_robust{robust}.txt', 'w') as f:
    f.write(f'# {{image_final}}.image\n')
    f.write(f'# nsigma: {nsigma}  fastnoise: False\n')
    f.write(f'# Beammajor(arcsec)  Beamminor(arcsec)  PA(deg)  Flux_in_mask(mJy)  Peak_int(mJy/beam)  rms(mJy/beam)  Peak_SNR\n')
    f.write(f'{{beammajor:.3f}}   {{beamminor:.3f}}   {{beampa:.2f}}   {{target_flux*1e3:.2f}}   {{peak_intensity*1e3:.2f}}   {{rms_mJy:.2e}}   {{SNR:.2f}}\n')

print(f"Finished imaging with robust {robust}")

print('\n##### EXTRACT UVTABLE #####\n')

tb.open("{msfile}")
data   = np.squeeze(tb.getcol("DATA"))
flag   = np.squeeze(tb.getcol("FLAG"))
uvw    = tb.getcol("UVW")
weight = tb.getcol("WEIGHT")
spwid  = tb.getcol("DATA_DESC_ID")
times  = tb.getcol("TIME")
tb.close()

tb.open("{msfile}" + "/SPECTRAL_WINDOW")
freqlist = np.squeeze(tb.getcol("CHAN_FREQ"))
if freqlist.ndim == 0:
    freqlist = np.array([freqlist])
tb.close()

good   = np.squeeze(np.any(flag, axis=0) == False)
data   = data[:, good]
weight = weight[:, good]
uvw    = uvw[:, good]
spwid  = spwid[good]
times  = times[good]

Re  = np.sum(data.real * weight, axis=0) / np.sum(weight, axis=0)
Im  = np.sum(data.imag * weight, axis=0) / np.sum(weight, axis=0)
vis = Re + 1j * Im
wgt = np.sum(weight, axis=0)

get_freq = lambda ispw: freqlist[ispw]
freqs = get_freq(spwid)

um, vm = uvw[0,:], uvw[1,:]
u = um * freqs / 2.9979e8
v = vm * freqs / 2.9979e8

os.makedirs("frank_fit", exist_ok=True)
if os.path.exists('frank_fit/uvtable_frank.npz'):
    os.remove('frank_fit/uvtable_frank.npz')
np.savez('frank_fit/uvtable_frank.npz', u=u, v=v, Vis=vis, Wgt=wgt)
print(f"Measurement set exported to frank_fit/uvtable_frank.npz")
