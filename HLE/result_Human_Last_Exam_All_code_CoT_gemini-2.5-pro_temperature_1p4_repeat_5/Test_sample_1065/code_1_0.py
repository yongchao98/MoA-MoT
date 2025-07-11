import math

# The center of mass coordinates are given by the equations:
# X_cm = (2/pi) * R
# Z_cm = (2/pi) * R
# We need to compute the numerical coefficient, which is 2/pi.

# The final equation for the coefficient is C = 2 / pi.
# To satisfy the instruction "output each number in the final equation!",
# we define the numbers from this equation as variables.
numerator = 2
denominator = math.pi

# Now, we calculate the coefficient using these numbers.
horizontal_coefficient = numerator / denominator
vertical_coefficient = numerator / denominator

# The problem asks for the raw numbers of the horizontal and vertical coordinates,
# separated by a comma. We print the calculated coefficients in this format.
print(f"{horizontal_coefficient},{vertical_coefficient}")