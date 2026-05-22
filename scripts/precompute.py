#!/usr/bin/env python3
"""
Precompute simulation data for the cerebral 0D model web visualization.

Runs the C simulator (./cbf) for 36 scenarios:
  2 conditions (NSR=0, AF=1) x 6 CoW variants (0-5) x 3 heart rates (60, 75, 90)

Extracts beats 100-109 (steady-state) from each run, normalizes time to start
at 0, and saves all scenarios to web/public/data/scenarios.json.
"""

import json
import os
import random
import subprocess
import sys

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CBF_BINARY = os.path.join(PROJECT_ROOT, "cbf")
INPUT_DIR = os.path.join(PROJECT_ROOT, "input")
OUTPUT_JSON = os.path.join(PROJECT_ROOT, "web", "public", "data", "scenarios.json")

CONDITIONS = [0, 1]  # 0 = NSR, 1 = AF
CONDITION_NAMES = {0: "nsr", 1: "af"}
COW_VARIANTS = list(range(6))
COW_NAMES = {
    0: "Complete (normal)",
    1: "Absent left PCoA",
    2: "Absent bilateral PCoA",
    3: "Absent left A1 (ACA)",
    4: "Absent left P1 (PCA)",
    5: "Absent right PCoA + left P1",
}
HEART_RATES = [60, 75, 90]

NUM_RANDOM_PARS = 95
NUM_BEATS = 6000
BEAT_START = 100  # first beat to extract (0-indexed)
BEAT_END = 109    # last beat to extract (inclusive, 0-indexed)

# The HR0 parameter is the 41st randomPar call (0-indexed: position 40)
HR_PAR_INDEX = 40
HR_BASE = 75.0  # base heart rate in the C code

COLUMNS = ["time", "P_a", "q_ml", "q_al", "q_pl", "q_mr", "q_ar", "q_pr"]

# ---------------------------------------------------------------------------
# Input file generation
# ---------------------------------------------------------------------------


def build_scenario_list():
    """Return a list of (run_index, condition, cow, hr) tuples."""
    scenarios = []
    run_index = 0
    for hr in HEART_RATES:
        for condition in CONDITIONS:
            for cow in COW_VARIANTS:
                scenarios.append((run_index, condition, cow, hr))
                run_index += 1
    return scenarios


def generate_input_files(scenarios):
    """
    Regenerate input/randomPars.dat, input/pinkNoise.dat, input/expNoise.dat
    with enough rows for all scenarios.

    - randomPars.dat: 95 tab-separated values of 1.0 per row, except column 40
      which is hr / 75.0 to set the desired heart rate.
    - pinkNoise.dat: 6000 tab-separated values of 0.0 per row.
    - expNoise.dat: 6000 tab-separated values from random.expovariate(1.0)
      with random.seed(42).
    """
    num_rows = len(scenarios)

    # --- randomPars.dat ---
    with open(os.path.join(INPUT_DIR, "randomPars.dat"), "w") as f:
        for run_index, _condition, _cow, hr in scenarios:
            values = [1.0] * NUM_RANDOM_PARS
            values[HR_PAR_INDEX] = hr / HR_BASE
            f.write("\t".join(str(v) for v in values) + "\n")

    # --- pinkNoise.dat ---
    with open(os.path.join(INPUT_DIR, "pinkNoise.dat"), "w") as f:
        zero_row = "\t".join(["0.0"] * NUM_BEATS)
        for _ in range(num_rows):
            f.write(zero_row + "\n")

    # --- expNoise.dat ---
    random.seed(42)
    with open(os.path.join(INPUT_DIR, "expNoise.dat"), "w") as f:
        for _ in range(num_rows):
            values = [random.expovariate(1.0) for _ in range(NUM_BEATS)]
            f.write("\t".join(str(v) for v in values) + "\n")

    print(f"Generated input files with {num_rows} rows each.")


# ---------------------------------------------------------------------------
# Simulation runner
# ---------------------------------------------------------------------------


def run_simulation(run_index, condition, cow, hr):
    """
    Run ./cbf run_index is_af cow_var hr_0 and return the output filenames.
    """
    cond_label = "AF" if condition == 1 else "NSR"
    key = f"{CONDITION_NAMES[condition]}_cow{cow}_hr{hr}"
    print(f"  Running scenario {key} (run_index={run_index})...", end="", flush=True)

    cmd = [CBF_BINARY, str(run_index), str(condition), str(cow), str(hr)]
    result = subprocess.run(cmd, cwd=PROJECT_ROOT, capture_output=True, text=True)
    if result.returncode != 0:
        print(f" FAILED")
        print(f"    stderr: {result.stderr.strip()}")
        sys.exit(1)

    state_file = os.path.join(
        PROJECT_ROOT, f"statesOutput.{run_index:05d}.{cond_label}.dat"
    )
    edt_file = os.path.join(
        PROJECT_ROOT, f"endDiastolicTime.{run_index:05d}.{cond_label}.dat"
    )
    param_file = os.path.join(PROJECT_ROOT, f"parameters.{run_index:05d}.dat")

    print(f" done.")
    return state_file, edt_file, param_file


# ---------------------------------------------------------------------------
# Data extraction
# ---------------------------------------------------------------------------


def extract_beats(state_file, edt_file, beat_start=BEAT_START, beat_end=BEAT_END):
    """
    Extract timesteps for beats beat_start through beat_end from the
    statesOutput file using the endDiastolicTime file.

    The end-diastolic time file has one time per line (one per beat).
    Beat N starts at edt[N-1] (end of beat N-1) and ends at edt[N].
    Beat 0 starts at time 0.

    Returns a list of rows, each row being [time, P_a, q_ml, q_al, q_pl, q_mr, q_ar, q_pr].
    Time is normalized to start at 0.
    """
    # Read end-diastolic times
    with open(edt_file, "r") as f:
        edt_times = [float(line.strip()) for line in f if line.strip()]

    # The start time of beat N:
    #   beat 0 starts at t=0
    #   beat N starts at edt_times[N-1] (the end-diastolic time of the previous beat)
    if beat_start == 0:
        t_start = 0.0
    else:
        t_start = edt_times[beat_start - 1]

    t_end = edt_times[beat_end]

    # Read state data and filter to time window
    rows = []
    with open(state_file, "r") as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) < 8:
                continue
            t = float(parts[0])
            if t < t_start:
                continue
            if t > t_end:
                break
            rows.append([float(x) for x in parts[:8]])

    if not rows:
        print(f"    WARNING: No data extracted for beats {beat_start}-{beat_end}")
        return []

    # Normalize time to start at 0
    t_offset = rows[0][0]
    for row in rows:
        row[0] -= t_offset

    # Round values: time to 4 decimals, others to 6 decimals
    for row in rows:
        row[0] = round(row[0], 4)
        for i in range(1, len(row)):
            row[i] = round(row[i], 6)

    return rows


# ---------------------------------------------------------------------------
# Cleanup
# ---------------------------------------------------------------------------


def cleanup_files(*files):
    """Remove output files produced by the simulator."""
    for f in files:
        if os.path.exists(f):
            os.remove(f)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------


def main():
    # Ensure output directory exists
    os.makedirs(os.path.dirname(OUTPUT_JSON), exist_ok=True)

    # Build scenario list
    scenarios = build_scenario_list()
    print(f"Will run {len(scenarios)} scenarios.")

    # Generate input files with enough rows
    generate_input_files(scenarios)

    # Run all simulations and extract data
    results = {}
    for run_index, condition, cow, hr in scenarios:
        key = f"{CONDITION_NAMES[condition]}_cow{cow}_hr{hr}"

        state_file, edt_file, param_file = run_simulation(
            run_index, condition, cow, hr
        )

        # Extract steady-state beats
        data = extract_beats(state_file, edt_file)
        print(f"    Extracted {len(data)} timesteps for beats {BEAT_START}-{BEAT_END}")

        results[key] = {
            "condition": CONDITION_NAMES[condition],
            "cow": cow,
            "hr": hr,
            "data": data,
        }

        # Clean up .dat output files
        cleanup_files(state_file, edt_file, param_file)

    # Build final JSON structure
    output = {
        "columns": COLUMNS,
        "cow_names": COW_NAMES,
        "scenarios": results,
    }

    with open(OUTPUT_JSON, "w") as f:
        json.dump(output, f, separators=(",", ":"))

    file_size_mb = os.path.getsize(OUTPUT_JSON) / (1024 * 1024)
    print(f"\nWrote {OUTPUT_JSON}")
    print(f"  File size: {file_size_mb:.2f} MB")
    print(f"  Scenarios: {len(results)}")
    print("Done.")


if __name__ == "__main__":
    main()
