import numpy as np
import os

print('\n##### RUNNING SIMOBSERVE #####\n')

os.chdir("output/ring10")  # Move into output folder

simobserve(
       project="gaussian_ring_B6_C-7_t1800s",
       skymodel="model_image.fits",
       indirection="",
       incell="",
       inbright="",
       incenter="230000000000.0Hz",
       inwidth="7500000000.0Hz",
       antennalist="../../alma-configuration-files/alma.cycle11.7.cfg",
       totaltime="1800s",
       user_pwv=1.5,
       overwrite=True
)


print(' ')
print('##### RUNNING TCLEAN #####')
print(' ')

# Define mask geometry

mask = f"ellipse[[16h00m00.000000s,-40d00m00.00000s], [0.656arcsec, 0.464arcsec], 30.0deg]"
noise_annulus = f"annulus[[16h00m00.000000s, -40d00m00.00000s],['1.094arcsec', '1.531arcsec']]"

# Estimate cell size based on resolution
#cell_arcsec = 0.015  # arcsec
imsize = 512         # adjusted for disk coverage

# Output paths
os.makedirs('tclean', exist_ok=True)
image_base = "tclean/ring10"

# First pass: low-niter clean to estimate RMS
image_first = f"{image_base}_firstclean"
for ext in [".image", ".mask", ".model", ".pb", ".psf", ".residual", ".sumwt", "weight"]:
    os.system(f"rm -rf {image_first}{ext}")

tclean(
    vis="gaussian_ring_B6_C-7_t1800s/gaussian_ring_B6_C-7_t1800s.alma.cycle11.7.noisy.ms",
    imagename=image_first,
    datacolumn='data',
    specmode='mfs',
    deconvolver='multiscale',
    scales=[0, 8, 15, 30, 80],
    gridder='standard',
    weighting='briggs',
    robust=0.5,
    gain=0.1,
    threshold='0.0Jy',
    cycleniter=300,
    cyclefactor=1,
    imsize=imsize,
    cell=f"0.015arcsec",
    niter=150,
    mask=mask,
    interactive=False,
    savemodel='none',
    nterms=1
)

rms_mJy = imstat(imagename=image_first + '.image', region=noise_annulus)['rms'][0] * 1e3
print(f"Estimated RMS: {rms_mJy:.3f} mJy")


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


SetWeights("gaussian_ring_B6_C-7_t1800s/gaussian_ring_B6_C-7_t1800s.alma.cycle11.7.noisy.ms", rms_mJy*1e-3)





# Second pass: full clean with threshold based on RMS
threshold_val = 2.0 * rms_mJy
threshold_str = f"{threshold_val:.2f}mJy"

image_final = f"{image_base}_robust0.5_clean"
for ext in [".image", ".mask", ".model", ".pb", ".psf", ".residual", ".sumwt", "weight"]:
    os.system(f"rm -rf {image_final}{ext}")

tclean(
    vis="gaussian_ring_B6_C-7_t1800s/gaussian_ring_B6_C-7_t1800s.alma.cycle11.7.noisy.ms",
    imagename=image_final,
    datacolumn='data',
    specmode='mfs',
    deconvolver='multiscale',
    scales=[0, 8, 15, 30],
    gridder='standard',
    weighting='briggs',
    robust=0.5,
    gain=0.1,
    threshold=threshold_str,
    cycleniter=300,
    cyclefactor=1,
    imsize=imsize,
    cell=f"0.015arcsec",
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
    os.system(f'rm -rf {image_first}{ext}')

# Save summary to text
with open(f'Info_image_ring10_robust0.5.txt', 'w') as Info_txt:
    Info_txt.write(f'# {image_final}.image\n')
    Info_txt.write(f'# Beammajor(arcsec)    Beamminor(arcsec)    PA(deg)    Flux_in_mask(mJy)    Peak_int(mJy/beam)    rms(mJy/beam)    Peak_SNR\n')
    Info_txt.write(f'{beammajor:.3f}   {beamminor:.3f}   {beampa:.2f}   {target_flux*1e3:.2f}   {peak_intensity*1e3:.2f}   {rms*1e3:.2e}   {SNR:.2f}')

print(f"Finished the imaging of the continuum with robust 0.5")

