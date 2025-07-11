# Number of fluorochromes/channels in the experiment
num_channels = 5

# --- Control Calculation ---

# 1. Unstained Control: Always one sample with just beads.
unstained_controls = 1

# 2. Single-Stain Compensation Controls: One for each channel.
compensation_controls = num_channels

# 3. Fluorescence Minus One (FMO) Controls: One for each channel to set accurate gates.
fmo_controls = num_channels

# --- Total Calculation ---
total_controls = unstained_controls + compensation_controls + fmo_controls

# --- Output the results ---
print(f"For a {num_channels}-color flow cytometry sorting experiment, you should prepare the following essential controls:")
print("-" * 20)
print(f"Unstained Controls:      {unstained_controls}")
print(f"Compensation Controls:   {compensation_controls} (one for each of your {num_channels} fluorochromes)")
print(f"FMO Controls:            {fmo_controls} (one for each fluorochrome to set accurate gates)")
print("-" * 20)
# Final equation as requested, showing each number
print(f"Total Essential Controls = {unstained_controls} (Unstained) + {compensation_controls} (Compensation) + {fmo_controls} (FMO) = {total_controls}")
