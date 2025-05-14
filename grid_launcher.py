import os
import json
import subprocess
from shutil import copyfile

# Base config to modify
base_config_path = "configs/1st_test_ring.json"
with open(base_config_path) as f:
    base_config = json.load(f)

# Grid of r1 values (ring location)
r1_values = [0, 5, 10, 20, 30, 40, 50]  # AU

# Where to save the new config files
os.makedirs("configs/grid_runs", exist_ok=True)

for r1 in r1_values:
    # Create a copy of the config with updated r1
    config = base_config.copy()
    config["r1"] = r1

    # Config name and path
    name = f"ring{r1}"
    config_filename = f"configs/grid_runs/{name}.json"

    # Write new config file
    with open(config_filename, "w") as f:
        json.dump(config, f, indent=4)

    # Run the full pipeline
    print(f"\nRunning pipeline for {name}")
    try:
        subprocess.run(["python", "run_pipeline.py", config_filename], check=True)
    except subprocess.CalledProcessError as e:
        print(f"Error running pipeline for {name}: {e.cmd}")
        continue
