# Step 1: Define the known values and fractions from the problem.
curves_stayed_true = 90
fraction_lost_shape = 3/8
fraction_new_paths = 1/4

# Step 2: Calculate the total fraction of curves that were affected.
# These are assumed to be two distinct groups.
total_fraction_affected = fraction_lost_shape + fraction_new_paths

# Step 3: Calculate the fraction of curves that stayed true.
# This is 1 minus the fraction of affected curves.
fraction_stayed_true = 1 - total_fraction_affected

# Step 4: Calculate the original total number of curves (T).
# T = (number of true curves) / (fraction of true curves)
original_total_curves = curves_stayed_true / fraction_stayed_true

# Step 5: Calculate the number of curves in each group for the final pattern.
# The number of unaffected curves is given.
unaffected_curves = curves_stayed_true

# The number of curves that lost their shape and must be redrawn.
redrawn_lost_shape_curves = original_total_curves * fraction_lost_shape

# The number of curves that found new paths and must be redrawn.
redrawn_new_path_curves = original_total_curves * fraction_new_paths

# Step 6: The total number of curves in the restored pattern is the sum of
# the unaffected curves and all the redrawn curves.
final_total_curves = unaffected_curves + redrawn_lost_shape_curves + redrawn_new_path_curves

# Step 7: Print the final calculation, showing how the total is composed.
# Each component must be an integer.
print("The final pattern is composed of the unaffected curves plus all the newly placed curves.")
print("The calculation for the total number of curves is:")
print(f"{int(unaffected_curves)} (unaffected) + {int(redrawn_lost_shape_curves)} (restored) + {int(redrawn_new_path_curves)} (rerouted) = {int(final_total_curves)}")
