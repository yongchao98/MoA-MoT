# Step 1: Establish the initial numbers from the problem description.
# The text states the original pattern had 360 curves.
# The poem's scenario is based on 90 curves remaining intact.
total_original_curves = 360
curves_remaining = 90

print(f"The original Rangoli pattern had {total_original_curves} curves.")
print(f"After the disruption, {curves_remaining} curves remained in their original form.")
print("-" * 20)

# Step 2: Calculate the total number of curves that were disturbed.
# These are the curves the master needs to redraw.
disturbed_curves = total_original_curves - curves_remaining

print("To find the number of curves the master must draw, we subtract the remaining curves from the original total.")
print(f"Calculation: {total_original_curves} (original) - {curves_remaining} (remaining) = {disturbed_curves} (to be drawn)")
print("-" * 20)

# Step 3: Verify the composition of the new curves based on the poem's fractions.
# The fractions apply to the total number of disturbed curves.
parabolic_curves = int((1/5) * disturbed_curves)
elliptical_curves = int((2/9) * disturbed_curves)
circular_curves = disturbed_curves - parabolic_curves - elliptical_curves

print("The poem describes the transformation of these disturbed curves:")
print(f"- Parabolic curves: (1/5) of {disturbed_curves} = {parabolic_curves}")
print(f"- Elliptical curves: (2/9) of {disturbed_curves} = {elliptical_curves}")
print(f"- Circular curves (the rest): {disturbed_curves} - {parabolic_curves} - {elliptical_curves} = {circular_curves}")
print("-" * 20)

# Step 4: State the final answer.
# The total number of curves to be drawn is the sum of these new types, which equals the total number of disturbed curves.
# The "golden ratio" affects the geometric properties (length) of the new curves, not the quantity required to complete the pattern.
total_curves_to_draw = disturbed_curves

print("The total number of curves the master must draw is the sum of all the new parabolic, elliptical, and circular curves.")
print(f"Final Equation: {parabolic_curves} + {elliptical_curves} + {circular_curves} = {total_curves_to_draw}")
print("-" * 20)
print(f"Final Answer: The master must draw a total of {total_curves_to_draw} curves to restore the pattern.")

# The final answer in the required format
print("<<<270>>>")