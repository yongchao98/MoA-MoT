# 1. Define the initial state based on the most consistent information.
initial_total_curves = 360
fraction_disturbed = 3/8

# 2. Calculate the number of curves the master needs to redraw.
# These are the curves that "lost their original form".
curves_to_redraw = int(initial_total_curves * fraction_disturbed)

# 3. Calculate the number of each new type of curve.
# The number of curves to redraw (135) is a multiple of lcm(5, 9)=45,
# which validates this approach.
parabolic_curves = int(curves_to_redraw * (1/5))
elliptical_curves = int(curves_to_redraw * (2/9))
circular_curves = curves_to_redraw - parabolic_curves - elliptical_curves

# 4. The total number of curves the master must draw is the sum of these new curves.
total_curves_to_draw = parabolic_curves + elliptical_curves + circular_curves

# Print the final breakdown and the total in an equation format.
# The question asks for the total number of curves the master must place.
print("To restore the pattern, the master must draw the new curves:")
print(f"Parabolic Curves: {parabolic_curves}")
print(f"Elliptical Curves: {elliptical_curves}")
print(f"Circular Curves: {circular_curves}")
print("\nThe final equation for the total curves to be drawn is:")
print(f"{parabolic_curves} + {elliptical_curves} + {circular_curves} = {total_curves_to_draw}")

# Final Answer
# print(f"<<<{total_curves_to_draw}>>>")