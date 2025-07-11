import math

# Step 1: Define the known values from the problem description and poem.
# The number of curves that remained undisturbed.
undisturbed_curves = 90

# The fractions of the total curves that were affected.
fraction_lost_shape = 3/8
fraction_found_new_paths = 1/4

# Step 2: Calculate the total fraction of curves that were affected by adding the two groups.
# This represents all curves that need to be replaced.
total_fraction_affected = fraction_lost_shape + fraction_found_new_paths

# Step 3: Calculate the fraction of curves that were undisturbed.
# This is the remainder of the total.
fraction_undisturbed = 1 - total_fraction_affected

# Step 4: Use the number of undisturbed curves to find the original total.
# The equation is: original_total * fraction_undisturbed = undisturbed_curves
# Therefore, original_total = undisturbed_curves / fraction_undisturbed
original_total_curves = undisturbed_curves / fraction_undisturbed

# Step 5: The question asks for the total number of curves in the restored pattern.
# A restored pattern must have the same number of curves as the original.
# The information about the golden ratio, symmetry, and specific curve types
# describes the characteristics of the new curves, not a change in their total quantity.
# The final answer is the original total number of curves.

# We will now print the components of the key calculation: Total = 90 / (3/8)
numerator = undisturbed_curves
denominator_part1 = 3
denominator_part2 = 8
final_answer = int(original_total_curves)

print("To find the total number of curves, we use the fact that 90 curves remained undisturbed.")
print(f"The fraction of affected curves is 3/8 + 1/4 = 5/8.")
print(f"Therefore, the fraction of undisturbed curves is 1 - 5/8 = 3/8.")
print(f"The final calculation is based on the equation: Total Curves * (3/8) = 90")
print(f"Solving for the total gives: Total = {numerator} / ({denominator_part1}/{denominator_part2})")
print(f"This simplifies to the final equation: ({numerator} * {denominator_part2}) / {denominator_part1} = {final_answer}")
print(f"\nThe total number of curves the master must have to restore the pattern is {final_answer}.")

# Verification of the breakdown of new curves (optional sanity check)
# This part confirms our interpretation is consistent.
# Total affected curves = total_fraction_affected * final_answer = (5/8) * 240 = 150
# "Lost shape" group (assumed to be the "disturbed" group for fraction breakdown) = fraction_lost_shape * final_answer = (3/8) * 240 = 90
# Parabolic curves = (1/5) * 90 = 18
# Elliptical curves = (2/9) * 90 = 20
# Circular curves = 90 - 18 - 20 = 52
# Since 18, 20, and 52 are all integers, our initial calculation of 240 is consistent with all constraints.

<<<240>>>