import math

# Step 1: Define the key numerical values from the problem description.
original_volume = 750  # The original sauce volume in ml
baseline_volume = 180  # The fixed baseline portion in ml
viscosity_increase_fraction = 3/7  # The fractional increase in viscosity

# Step 2: Calculate the original adjustable portion of the sauce.
# This is the part of the volume that can be reduced due to increased efficiency.
adjustable_volume = original_volume - baseline_volume

# Step 3: Calculate the new, reduced adjustable volume.
# The volume reduction is modeled with an exponential decay function, which corresponds to the
# required "logarithmic relationship" with the viscosity increase.
reduction_factor = math.exp(-viscosity_increase_fraction)
new_adjustable_volume = adjustable_volume * reduction_factor

# Step 4: Calculate the final total volume for the new sauce.
# This is the sum of the fixed baseline volume and the new reduced adjustable volume.
final_volume = baseline_volume + new_adjustable_volume

# Step 5: Print the final answer, showing the full equation as requested.
# The numbers from the problem are explicitly included in the output string.
print("To find the required volume of the new sauce, we model the change based on the provided baseline and viscosity increase.")
print("The final volume is the sum of the fixed baseline and a reduced adjustable portion.\n")
print(f"Final Volume = Baseline Volume + (Original Volume - Baseline Volume) * e^(-Viscosity Increase)")
print(f"Final Volume = {baseline_volume} ml + ({original_volume} ml - {baseline_volume} ml) * e^(-3/7)")
print(f"Final Volume = {baseline_volume} ml + {adjustable_volume} ml * {reduction_factor:.4f}")
print(f"Final Volume = {baseline_volume} ml + {new_adjustable_volume:.1f} ml")
print(f"Total New Sauce Required = {final_volume:.1f} ml")
