# Step 1: Define the initial values from the problem description.
total_original_curves = 360
curves_unchanged = 90

print(f"The original Rangoli pattern had {total_original_curves} curves.")
print(f"After the disruption, {curves_unchanged} curves remained unchanged.")
print("-" * 30)

# Step 2: Calculate the number of curves that are replaced.
# The poem states "one-fourth find new pathways," which are the curves being replaced.
fraction_finding_new_paths = 1/4
new_path_curves = int(total_original_curves * fraction_finding_new_paths)

print(f"The number of curves that found new pathways is 1/4 of the total: {new_path_curves}")

# Step 3: Calculate the breakdown of the new curves.
# "Of those that left their plotted way" refers to these {new_path_curves} curves.
# - One-fifth become parabolic
# - Two-ninths become elliptical
# - The rest become circular

parabolic_fraction = 1/5
elliptical_fraction = 2/9

parabolic_curves = int(new_path_curves * parabolic_fraction)
elliptical_curves = int(new_path_curves * elliptical_fraction)
circular_curves = new_path_curves - parabolic_curves - elliptical_curves

print("\nThese new curves are composed of:")
print(f" - Parabolic curves (1/5 of {new_path_curves}): {parabolic_curves}")
print(f" - Elliptical curves (2/9 of {new_path_curves}): {elliptical_curves}")
print(f" - Circular curves (the remainder): {circular_curves}")
print("-" * 30)

# Step 4: Calculate the total number of curves in the restored pattern.
# The restored pattern consists of the unchanged curves plus the newly drawn curves.
# The golden ratio affects curve length, not the count, so it does not alter this calculation.
total_final_curves = curves_unchanged + new_path_curves

print("To restore the pattern, the master must have the original unchanged curves plus all the new curves.")
print("Total Curves in Restored Pattern = Unchanged Curves + New Curves")
print("\nFinal Equation:")
print(f"{curves_unchanged} + {new_path_curves} = {total_final_curves}")
