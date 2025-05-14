import numpy as np
import os

print('\n##### RUNNING SIMOBSERVE #####\n')

os.chdir("{output_dir}")  # Move into output folder

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


print(' ')
print('##### RUNNING TCLEAN #####')
print(' ')

# Define mask geometry

mask = f"ellipse[[{mask_ra},{mask_dec}], [{mask_semimajor:.3f}arcsec, {mask_semiminor:.3f}arcsec], {PA:.1f}deg]"
noise_annulus = f"annulus[[{mask_ra}, {mask_dec}],['{noise_inner:.3f}arcsec', '{noise_outer:.3f}arcsec']]"

# Estimate cell size based on resolution
#cell_arcsec = {cellsize}  # arcsec
imsize = {imsize}         # adjusted for disk coverage

# Output paths
os.makedirs('tclean', exist_ok=True)
image_base = "tclean/{base_name}"

# First pass: low-niter clean to estimate RMS
image_first = f"{{image_base}}_firstclean"
for ext in [".image", ".mask", ".model", ".pb", ".psf", ".residual", ".sumwt", "weight"]:
    os.system(f"rm -rf {{image_first}}{{ext}}")

tclean(
    vis="{msfile}",
    imagename=image_first,
    datacolumn='data',
    specmode='mfs',
    deconvolver='multiscale',
    scales=[0, 8, 15, 30, 80],
    gridder='standard',
    weighting='briggs',
    robust={robust},
    gain=0.1,
    threshold='0.0Jy',
    cycleniter=300,
    cyclefactor=1,
    imsize=imsize,
    cell=f"{cellsize:.3f}arcsec",
    niter=150,
    mask=mask,
    interactive=False,
    savemodel='none',
    nterms=1
)

rms_mJy = imstat(imagename=image_first + '.image', region=noise_annulus)['rms'][0] * 1e3
print(f"Estimated RMS: {{rms_mJy:.3f}} mJy")


# Setting the right weights
def SetWeights(MS_filename, rms):

     # rms is the rms noise in Jy
     tb.open(MS_filename, nomodify=False)
     weight = tb.getcol("WEIGHT")
     sigma = tb.getcol("SIGMA")

     nvis = weight.shape[1]
     weight[:]=1./2./nvis/rms**2
     sigma = 1./np.sqrt(weight)

     tb.putcol("WEIGHT", weight)
     tb.putcol("SIGMA", sigma)

     tb.flush()
     tb.close()


SetWeights("{msfile}", rms_mJy*1e-3)





# Second pass: full clean with threshold based on RMS
threshold_val = {threshold} * rms_mJy
threshold_str = f"{{threshold_val:.2f}}mJy"

image_final = f"{{image_base}}_robust{robust}_clean"
for ext in [".image", ".mask", ".model", ".pb", ".psf", ".residual", ".sumwt", "weight"]:
    os.system(f"rm -rf {{image_final}}{{ext}}")

tclean(
    vis="{msfile}",
    imagename=image_final,
    datacolumn='data',
    specmode='mfs',
    deconvolver='multiscale',
    scales=[0, 8, 15, 30],
    gridder='standard',
    weighting='briggs',
    robust={robust},
    gain=0.1,
    threshold=threshold_str,
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
exportfits(imagename=image_final + '.image',
           fitsimage=image_final + '.image.fits',
           overwrite=True, history=True)

# Measure beam and source properties
headerlist = imhead(image_final + '.image', mode='list')
beammajor = headerlist['beammajor']['value']
beamminor = headerlist['beamminor']['value']
beampa = headerlist['beampa']['value']
target_stats = imstat(imagename=image_final + '.image', region=mask)
target_flux = target_stats['flux'][0]
peak_intensity = target_stats['max'][0]
rms = imstat(imagename=image_final + '.image', region=noise_annulus)['rms'][0]
SNR = peak_intensity / rms

# Delete firstclean results
for ext in ['.image', '.mask', '.model', '.pb', '.psf', '.residual', '.sumwt', 'weight']:
    os.system(f'rm -rf {{image_first}}{{ext}}')

# Save summary to text
with open(f'Info_image_{base_name}_robust{robust}.txt', 'w') as Info_txt:
    Info_txt.write(f'# {{image_final}}.image\n')
    Info_txt.write(f'# Beammajor(arcsec)    Beamminor(arcsec)    PA(deg)    Flux_in_mask(mJy)    Peak_int(mJy/beam)    rms(mJy/beam)    Peak_SNR\n')
    Info_txt.write(f'{{beammajor:.3f}}   {{beamminor:.3f}}   {{beampa:.2f}}   {{target_flux*1e3:.2f}}   {{peak_intensity*1e3:.2f}}   {{rms*1e3:.2e}}   {{SNR:.2f}}')

print(f"Finished the imaging of the continuum with robust {robust}")

