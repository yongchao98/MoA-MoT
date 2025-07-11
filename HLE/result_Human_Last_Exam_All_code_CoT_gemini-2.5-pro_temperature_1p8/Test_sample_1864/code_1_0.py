# Number of fluorochromes specified in the experiment
num_colors = 5

# 1. Unstained control: A sample with just the beads to measure autofluorescence.
unstained_controls = 1

# 2. Single-Color Compensation Controls: One control for each fluorochrome.
# This is necessary to create the compensation matrix to correct for spectral overlap.
compensation_controls = num_colors

# 3. Fluorescence Minus One (FMO) Controls: One FMO for each channel.
# These are used to accurately set gates by showing the spread of all other
# fluorochromes into the channel of interest.
fmo_controls = num_colors

# Calculate the total number of essential controls
total_controls = unstained_controls + compensation_controls + fmo_controls

# Print the breakdown and the final equation
print("Calculating the number of essential technical controls for a 5-color sorting experiment:\n")
print(f"1. Unstained Control: This is a single sample for baseline fluorescence.")
print(f"   - Count: {unstained_controls}\n")
print(f"2. Single-Color Compensation Controls: One is needed for each of the {num_colors} colors.")
print(f"   - Count: {compensation_controls}\n")
print(f"3. Fluorescence Minus One (FMO) Controls: One is needed for each of the {num_colors} colors to set accurate gates.")
print(f"   - Count: {fmo_controls}\n")
print("------------------------------------------------------------------")
print("Total number of controls is the sum of these types.")
print(f"Final Equation: {unstained_controls} (Unstained) + {compensation_controls} (Compensation) + {fmo_controls} (FMO) = {total_controls}")
print("------------------------------------------------------------------")
print(f"\nYou should prepare {total_controls} essential technical control tubes for your experiment.")