import math

# Step 1: Define the diameters of the two spaces.
# The diameter of the interval [0, 1] is the distance between its endpoints.
diameter_interval = 1.0

# The diameter of the unit circle with the intrinsic metric is half its circumference.
# Circumference C = 2 * pi * r. With radius r=1, C = 2 * pi.
# Diameter = (2 * pi) / 2 = pi.
diameter_circle = math.pi

# Step 2: Apply the formula for the Gromov-Hausdorff distance.
# For these specific spaces, the distance is given by the formula:
# d_GH = (1/2) * |diameter_circle - diameter_interval|
# Since pi > 1, the absolute value is not necessary.
distance = (diameter_circle - diameter_interval) / 2

# Step 3: Print the final equation and the result.
# The user wants to see each number in the final equation.
pi_val = math.pi
numerator_val1 = pi_val
numerator_val2 = diameter_interval
denominator_val = 2.0

print(f"The Gromov-Hausdorff distance is calculated by the formula: (diam(S^1) - diam([0,1])) / 2")
print(f"Final Equation: ({numerator_val1} - {numerator_val2}) / {denominator_val}")
print(f"Result: {distance}")