import math

# The horizontal and vertical coordinates of the center of mass are (2R/π, 2R/π).
# The problem asks for the raw numerical coefficients of R in the final equation for the coordinates.
# The final equation for the horizontal coordinate is X_cm = (2/π) * R.
# The final equation for the vertical coordinate is Y_cm = (2/π) * R.
# We will calculate the value of the coefficient 2/π.

# The numbers in the final equation are 2 and pi.
numerator = 2
denominator = math.pi

# Calculate the numerical value of the coefficient.
coefficient = numerator / denominator

# Print the horizontal and vertical coordinates separated by a comma.
print(f"{coefficient},{coefficient}")