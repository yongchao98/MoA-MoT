import math

# The plan is outlined in the thinking steps above.
# The final coordinates for the center of mass are (X_cm, Y_cm).
# X_cm = 0
# Y_cm = R * (1 - pi / 4)
# The problem asks for the raw numbers, so we set the radius R to 1.

R = 1.0

# Define the numbers used in the final equations.
x_val = 0.0
pi_val = math.pi
one = 1.0
four = 4.0

# The horizontal coordinate is 0, as the string hangs vertically along the y-axis.
x_cm = x_val

# The vertical coordinate is calculated from the equation Y_cm = R * (1 - pi / 4).
# We print the components of the equation first as requested.
print(f"The equation for the vertical coordinate is Y_cm = R * ({one} - {pi_val} / {four})")
y_cm = R * (one - pi_val / four)

# Finally, print the raw numerical coordinates separated by a comma.
print(f"The calculated coordinates (x, y) are:")
print(f"{x_cm},{y_cm}")