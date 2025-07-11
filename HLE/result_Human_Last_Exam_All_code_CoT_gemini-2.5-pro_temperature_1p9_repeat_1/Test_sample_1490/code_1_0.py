import math

# Step 1: Define the known values from the problem description.
# The number of curves that were not disturbed.
curves_stay_true = 90
# The fraction of curves that lost their original form.
fraction_lost_form = 3/8
# The fraction of curves that found new pathways.
fraction_new_pathways = 1/4

# Step 2: Calculate the total fraction of disturbed curves and the fraction that remained true.
total_disturbed_fraction = fraction_lost_form + fraction_new_pathways
fraction_true = 1 - total_disturbed_fraction

# Step 3: Calculate the total number of curves in the original pattern.
# We know that `fraction_true` of the total corresponds to `curves_stay_true`.
# So, total_curves = curves_stay_true / fraction_true.
if fraction_true > 0:
    total_initial_curves = int(curves_stay_true / fraction_true)
else:
    total_initial_curves = 0

# Step 4: To restore the pattern, all disturbed curves must be redrawn.
# This means the total number of curves in the final pattern is the same as the initial total.
# We can calculate the number of curves in each group that was redrawn.
curves_redrawn_lost_form = int(fraction_lost_form * total_initial_curves)
curves_redrawn_new_pathways = int(fraction_new_pathways * total_initial_curves)

# The total number of curves in the final, restored pattern is the sum of the
# unchanged curves and all the redrawn curves.
total_final_curves = curves_stay_true + curves_redrawn_lost_form + curves_redrawn_new_pathways

# Step 5: Display the final equation and the answer.
print("To find the total number of curves in the restored pattern, we sum the parts:")
print(f"Unchanged Curves + Redrawn 'Lost Form' Curves + Redrawn 'New Pathway' Curves = Total Curves")
print("The final equation is:")
print(f"{curves_stay_true} + {curves_redrawn_lost_form} + {curves_redrawn_new_pathways} = {total_final_curves}")
print("\nThus, the total number of curves in the restored pattern is:")
print(total_final_curves)
