# DustForge: Simulating ALMA Dust Substructure

DustForge is a modular pipeline for simulating ALMA continuum observations of protoplanetary disks with embedded dust substructures.

## Features
- Define analytical disk intensity profiles (e.g., Gaussian rings)
- Convert radial profiles into 2D images with inclination and position angle
- Normalize to desired total flux and disk size
- Simulate realistic ALMA observations with `simobserve`
- Perform CLEAN imaging with `tclean`, extract SNR and beam stats
- Flexible JSON-based configuration
- Easy to scale across multiple parameter grids

##  Requirements

- Python 3.7+
- CASA 6+ (monolithic version)
- `astropy`, `numpy`, `matplotlib`

## Folder structure
```
DustForge/
├── build_model_image.py          # Generates 2D disk image from radial profile
├── simulate_observation.py       # Runs simobserve + tclean
├── CASA_scripts.py               # CASA commands template
├── ALMA_specs.py                 # Default ALMA band/config specs
├── run_pipeline.py               # Runs model + synthetic observation in one step
├── configs/                      # Input JSON config files
├── models/                       # Intensity profile functions
├── output/                       # All outputs go here
└── README.md
```

## Quickstart
```
# Run the full pipeline in one step
python run_pipeline.py configs/1st_test_ring.json
```
Output products will be saved under output/{config_name}/ and include:
- `model_image.fits`: input sky model
- `tclean/*.image.fits`: cleaned image
- `Info_image_*.txt`: beam, flux, RMS, SNR
