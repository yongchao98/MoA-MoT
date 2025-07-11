# 1. Define initial known values from the problem.
original_total_volume_ml = 750.0
baseline_portion_ml = 180.0
eggs = 6
original_displacement_per_egg_ml = 45.0
egg_size_increase_factor = 1.12  # 12% larger
viscosity_increase_fraction_numerator = 3
viscosity_increase_fraction_denominator = 7

# 2. Calculate the "adjustable" portion of the original sauce.
adjustable_portion_ml = original_total_volume_ml - baseline_portion_ml

# 3. Calculate the volume efficiency factor based on the viscosity change.
# New viscosity = (1 + 3/7) = 10/7 times the old viscosity.
# The required volume is inversely proportional to the viscosity.
volume_efficiency_factor = 1 / (1 + viscosity_increase_fraction_numerator / viscosity_increase_fraction_denominator)

# Calculate the new volume of the adjusted portion.
new_adjusted_volume_ml = adjustable_portion_ml * volume_efficiency_factor

# 4. Calculate the change in volume due to larger eggs.
original_total_displacement_ml = eggs * original_displacement_per_egg_ml
new_total_displacement_ml = original_total_displacement_ml * egg_size_increase_factor
extra_displacement_volume_ml = new_total_displacement_ml - original_total_displacement_ml

# 5. Calculate the final total sauce volume for the new recipe.
final_volume_ml = baseline_portion_ml + new_adjusted_volume_ml + extra_displacement_volume_ml

# 6. Print the final equation and the result, showing each component.
print("To maintain the critical ratios with the new ingredients, Chef Sarah must prepare a new volume of sauce.")
print("The calculation is as follows:")
print(f"Final Volume = Baseline Portion + Adjusted Sauce Portion + Extra Egg Displacement")
print(f"Final Volume = {baseline_portion_ml:.2f} ml + {new_adjusted_volume_ml:.2f} ml + {extra_displacement_volume_ml:.2f} ml")
print(f"Total New Sauce Volume = {final_volume_ml:.2f} ml")

# Final answer in the required format
# print(f"<<<{final_volume_ml:.1f}>>>")