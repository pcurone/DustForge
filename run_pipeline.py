import sys
import os
import subprocess

# --- Select which steps to run ---
RUN_MODEL = False
RUN_SIMULATION = False
RUN_ANALYSIS = True

if len(sys.argv) < 2:
    print("Usage: python run_pipeline.py <config_file.json>")
    sys.exit(1)

config = sys.argv[1]

try:
    if RUN_MODEL:
        subprocess.run(["python", "build_model_image.py", config], check=True)
    if RUN_SIMULATION:
        subprocess.run(["python", "simulate_observation.py", config], check=True)
    if RUN_ANALYSIS:
        subprocess.run(["python", "recovery_score.py", config], check=True)
except subprocess.CalledProcessError as e:
    print(f"\nPipeline stopped due to an error in: {e.cmd}")
    sys.exit(1)