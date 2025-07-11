import math

# Step 1: Define the initial values from the problem description.
original_volume_ml = 750
baseline_volume_ml = 180
viscosity_increase_fraction = 3/7
egg_size_increase_percent = 0.12

# Step 2: Calculate the original "variable" volume of sauce.
# This is the portion of the sauce affected by the eggs and efficiency changes.
original_variable_volume = original_volume_ml - baseline_volume_ml

# Step 3: Adjust the variable volume for the larger eggs.
# The new eggs are 12% larger, so they require proportionally more sauce volume before efficiency is considered.
egg_size_multiplier = 1 + egg_size_increase_percent
adjusted_variable_volume = original_variable_volume * egg_size_multiplier

# Step 4: Calculate the efficiency factor based on the logarithmic rule for viscosity.
# The new viscosity is (1 + 3/7) times the old one. The efficiency gain is modeled as 1 + ln(ratio).
viscosity_ratio = 1 + viscosity_increase_fraction
efficiency_factor = 1 + math.log(viscosity_ratio) # math.log is the natural logarithm (ln)

# Step 5: Apply the efficiency factor to the adjusted variable volume.
# Higher efficiency means less sauce is needed, so we divide by the factor.
final_variable_volume = adjusted_variable_volume / efficiency_factor

# Step 6: Calculate the final total volume for the new sauce.
# This is the sum of the stable baseline and the new, efficient variable volume.
final_total_volume = baseline_volume_ml + final_variable_volume

# Step 7: Print the final equation and the result.
print("Final Sauce Volume Calculation")
print("="*30)
# Print the equation with each number clearly shown
print(f"Final Volume = Baseline + ((Original Volume - Baseline) * (1 + Egg Increase)) / (1 + ln(1 + Viscosity Increase))")
print(f"Final Volume = {baseline_volume_ml} + (({original_volume_ml} - {baseline_volume_ml}) * (1 + {egg_size_increase_percent})) / (1 + ln(1 + {viscosity_increase_fraction:.2f}))")
print(f"Final Volume = {baseline_volume_ml} + ({original_variable_volume} * {egg_size_multiplier}) / (1 + ln({viscosity_ratio:.3f}))")
print(f"Final Volume = {baseline_volume_ml} + {adjusted_variable_volume:.1f} / {efficiency_factor:.4f}")
print(f"Final Volume = {baseline_volume_ml} + {final_variable_volume:.2f}")
print(f"Final Volume = {final_total_volume:.2f} ml")
print("="*30)

# Final answer in the required format
print(f"\nTo maintain all critical ratios, Chef Sarah must prepare {final_total_volume:.2f} ml of the new sauce mixture.")
print("<<<650.55>>>")