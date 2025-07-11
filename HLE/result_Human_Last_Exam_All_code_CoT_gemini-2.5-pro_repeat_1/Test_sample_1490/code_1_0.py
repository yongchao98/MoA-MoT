# Rangoli Pattern Restoration Calculation

# Step 1: Define the initial conditions from the problem description.
initial_total_curves = 360
curves_staying_true = 90

# Step 2: Calculate the number of curves that were disturbed and must be redrawn.
# The problem provides conflicting information (3/8 of 360 disturbed vs. 90 staying true).
# We use the explicit number '90' as it is stated as a direct fact ("Given that: 90 curves...").
# This implies that the 'three-eighths' fraction is a poetic distractor.
curves_to_redraw = initial_total_curves - curves_staying_true

# Step 3: Calculate the breakdown of the newly formed curves.
# The fractions for new curve types apply to the set of curves being redrawn.
new_parabolic = int((1/5) * curves_to_redraw)
new_elliptical = int((2/9) * curves_to_redraw)
new_circular = curves_to_redraw - (new_parabolic + new_elliptical)

# Step 4: Calculate the total number of curves in the final, restored pattern.
# The restored pattern consists of the curves that were never changed plus all the new curves.
final_total_curves = curves_staying_true + curves_to_redraw

# Step 5: Print the step-by-step logic and the final answer.
# The output will show each number used in the final calculation.

print("--- Rangoli Restoration Analysis ---\n")
print(f"1. Initial State: The pattern starts with {initial_total_curves} curves.")
print(f"2. Unchanged Curves: It is given that {curves_staying_true} curves remained true to their original form.")
print(f"3. Curves to Restore: The number of curves the master must redraw is the difference:")
print(f"   {initial_total_curves} (Initial) - {curves_staying_true} (Unchanged) = {curves_to_redraw} curves.\n")

print(f"4. Breakdown of the {curves_to_redraw} New Curves:")
print(f"   - Parabolic: (1/5) of {curves_to_redraw} = {new_parabolic}")
print(f"   - Elliptical: (2/9) of {curves_to_redraw} = {new_elliptical}")
print(f"   - Circular: {curves_to_redraw} - ({new_parabolic} + {new_elliptical}) = {new_circular}\n")

print("--- Final Calculation ---")
print("The total number of curves in the fully restored pattern is the sum of the unchanged curves and all the redrawn curves.")
print("The final equation is:\n")
print(f"Total Curves = (Curves Staying True) + (Total Redrawn Curves)")
print(f"Total Curves = {curves_staying_true} + {curves_to_redraw} = {final_total_curves}")
print("\nOr, showing the full breakdown:")
print(f"Total Curves = {curves_staying_true} (Original) + {new_parabolic} (Parabolic) + {new_elliptical} (Elliptical) + {new_circular} (Circular) = {final_total_curves}")
print("\nTherefore, the total number of curves the master must draw to completely restore the pattern is 360.")

print("\n<<<360>>>")