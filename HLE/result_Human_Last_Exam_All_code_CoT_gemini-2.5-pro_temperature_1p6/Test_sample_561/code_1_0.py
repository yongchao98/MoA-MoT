import math

# Define the parameters of the self-affine fractal based on the problem description.
# N is the number of self-similar pieces (the number of black keys in an octave).
N = 5

# r_x is the scaling factor in the x-direction (width).
# Total keyboard width is 3.
# A white key is 3/7 units wide. A black key is half that, so 3/14 units wide.
# r_x = (width of black key) / (width of keyboard)
r_x = (3 / 14) / 3

# r_y is the scaling factor in the y-direction (height).
# Total keyboard height is 1. Black key height is 9/14.
# r_y = (height of black key) / (height of keyboard)
r_y = 9 / 14

# The formula for the box-counting dimension D for a self-affine fractal
# where r_x < r_y is D = log(N * r_y / r_x) / log(1 / r_x).

# Calculate the terms for the formula.
numerator_val = N * r_y / r_x
denominator_val = 1 / r_x

# Calculate the final dimension.
dimension = math.log(numerator_val) / math.log(denominator_val)

# Print the final equation with all the numbers.
print("The Minkowski–Bouligand dimension (D) is calculated using the formula for self-affine fractals.")
print(f"The equation is: D = log(N * r_y / r_x) / log(1 / r_x)")
print("\nPlugging in the values derived from the problem:")
print(f"N = {N}")
print(f"r_y / r_x = ({r_y}) / ({r_x}) = {r_y/r_x}")
print(f"1 / r_x = 1 / ({r_x}) = {1/r_x}")
print("\nThis simplifies the equation to:")
print(f"D = log({N} * {r_y/r_x}) / log({1/r_x})")
print(f"D = log({numerator_val}) / log({denominator_val})")
print(f"\nThe calculated dimension is: {dimension}")