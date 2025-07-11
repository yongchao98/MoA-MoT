import math

# --- Given values from the problem description ---
original_volume_ml = 750
baseline_volume_ml = 180
viscosity_increase_num = 3
viscosity_increase_den = 7
num_eggs = 6
original_egg_displacement_cm3 = 45
egg_size_increase_percent = 0.12

# --- Step 1: Calculate the volume affected by the viscosity change ---
# This is the portion of the sauce that is not part of the stable baseline.
affected_volume = original_volume_ml - baseline_volume_ml

# --- Step 2: Calculate the volume reduction due to increased viscosity ---
# The viscosity increases by 3/7, so the factor is 1 + 3/7 = 10/7.
# The problem states the relationship is logarithmic.
viscosity_factor = 1 + viscosity_increase_num / viscosity_increase_den
logarithmic_reduction = affected_volume * math.log(viscosity_factor)

# --- Step 3: Calculate the extra volume needed for the larger eggs ---
# The eggs are 12% larger, so we calculate the additional volume displaced by all 6 eggs.
extra_volume_for_eggs = num_eggs * original_egg_displacement_cm3 * egg_size_increase_percent

# --- Step 4: Calculate the final required sauce volume ---
# Start with the original volume, subtract the viscosity savings, and add the volume for larger eggs.
new_sauce_volume = original_volume_ml - logarithmic_reduction + extra_volume_for_eggs

# --- Step 5: Output the final equation and the result ---
# The prompt requires printing the equation with all the numbers.
print("The final calculation is based on adjusting the original volume:")
print("New Volume = Original Volume - Viscosity Reduction + Egg Volume Increase\n")

print("Equation:")
print(f"New Volume = {original_volume_ml}ml - ({original_volume_ml}ml - {baseline_volume_ml}ml) * log(1 + {viscosity_increase_num}/{viscosity_increase_den}) + ({num_eggs} * {original_egg_displacement_cm3}cm³ * {egg_size_increase_percent})")
print(f"New Volume = {original_volume_ml} - ({affected_volume}) * log({viscosity_factor:.3f}) + ({extra_volume_for_eggs})")
print(f"New Volume = {original_volume_ml} - {logarithmic_reduction:.3f} + {extra_volume_for_eggs}")
print(f"\nFinal calculated sauce volume: {new_sauce_volume:.1f} ml")

print("\n<<<579.1>>>")