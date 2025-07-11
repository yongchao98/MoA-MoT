# The initial number of curves in the Rangoli pattern.
initial_curves = 360

# The fraction for the curves that were disturbed. These are ultimately replaced,
# so they don't change the base count of the restored pattern.
# disturbed_curves = int((3/8) * initial_curves)

# The fraction for the new curves that were added to the pattern.
# These are described as "one-fourth found new paths to trace".
added_curves_fraction_numerator = 1
added_curves_fraction_denominator = 4

# Calculate the number of new curves added to the pattern.
added_curves = int(initial_curves * (added_curves_fraction_numerator / added_curves_fraction_denominator))

# The final total number of curves is the original number plus the newly added curves.
final_total_curves = initial_curves + added_curves

# Print the final equation showing each number, as requested.
# The final pattern consists of the original 360 curves (after restoration)
# plus the 90 new curves that were added.
print("The final pattern is composed of the restored original curves plus the newly added curves.")
print("The equation for the final total number of curves is:")
print(f"{initial_curves} + {added_curves} = {final_total_curves}")
