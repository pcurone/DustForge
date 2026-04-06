"""
run_pipeline.py
===============
Runs the full DustForge pipeline for a given JSON config:
    1. Build model image (build_model_image.py)
    2. Simulate observation — dispatches automatically:
         - simulate_observation_multiconfig.py  if "configurations" list is present
         - simulate_observation.py              for single-config runs
    3. Recovery analysis (recovery_score.py)

Usage
-----
    python run_pipeline.py configs/my_config.json
"""

import sys
import os
import json
import subprocess

# --- Select which steps to run ---
RUN_MODEL      = True
RUN_SIMULATION = True
RUN_ANALYSIS   = True

if len(sys.argv) < 2:
    print("Usage: python run_pipeline.py <config_file.json>")
    sys.exit(1)

config_path = sys.argv[1]

# Peek at the config to decide which simulation driver to use
with open(config_path) as f:
    config = json.load(f)

is_multiconfig = "configurations" in config

if is_multiconfig:
    sim_script = "simulate_observation_multiconfig.py"
    print(f"[run_pipeline] Multi-config mode detected "
          f"({len(config['configurations'])} configurations)")
else:
    sim_script = "simulate_observation.py"
    print(f"[run_pipeline] Single-config mode "
          f"({config.get('ALMA_configuration', '?')})")

try:
    if RUN_MODEL:
        subprocess.run(["python", "build_model_image.py", config_path], check=True)
    if RUN_SIMULATION:
        subprocess.run(["python", sim_script, config_path], check=True)
    if RUN_ANALYSIS:
        subprocess.run(["python", "recovery_score.py", config_path], check=True)
except subprocess.CalledProcessError as e:
    print(f"\nPipeline stopped due to an error in: {e.cmd}")
    sys.exit(1)
