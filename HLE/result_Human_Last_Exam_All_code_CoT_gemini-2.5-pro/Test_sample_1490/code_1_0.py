# The master needs to restore the pattern. The final pattern will be composed of:
# 1. The curves that remained unchanged.
# 2. The new curves drawn to replace those that "lost their original form".
# 3. The new curves drawn to replace those that "found new pathways".

# We are given the number of curves that stayed true.
curves_stayed_true = 90

# Let T be the original total number of curves. The problem states that the total (T)
# is the sum of those that stayed true (90), those that lost their form (3/8 * T),
# and those that found new paths (1/4 * T).
# The equation is: 90 + (3/8 * T) + (1/4 * T) = T
# This simplifies to: 90 = (3/8) * T
# Solving for T gives the original total number of curves.
original_total_curves = int(90 * (8 / 3))

# Now we can calculate the number of curves in each group that was disturbed and needs to be redrawn.
curves_lost_form = int((3/8) * original_total_curves)
curves_found_new_paths = int((1/4) * original_total_curves)

# The total number of curves in the restored pattern is the sum of the curves
# that stayed plus the new curves that were drawn to replace the disturbed ones.
final_total_curves = curves_stayed_true + curves_lost_form + curves_found_new_paths

# Print the final equation showing the composition of the restored pattern.
print("To restore the sacred symmetry, the master must ensure the final pattern contains:")
print(f"The {curves_stayed_true} curves that stayed true,")
print(f"plus the {curves_lost_form} curves redrawn where shapes were lost,")
print(f"plus the {curves_found_new_paths} curves redrawn along new paths.")
print("\nThe final equation for the total number of curves is:")
print(f"{curves_stayed_true} + {curves_lost_form} + {curves_found_new_paths} = {final_total_curves}")
print(f"\nThus, the master must draw a total of {final_total_curves} curves to restore the pattern.")
