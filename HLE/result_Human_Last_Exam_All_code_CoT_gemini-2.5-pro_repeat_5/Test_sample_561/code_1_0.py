import math

# Step 1: Define the parameters of the fractal construction.

# N is the number of self-similar pieces at each step of the construction.
# In one octave, there are 5 black keys.
N = 5

# r_x is the scaling factor in the horizontal (x) direction.
# The original keyboard width is 3 units.
# The width of a white key is 3/7 units.
# The width of a black key is half of a white key's width: (1/2) * (3/7) = 3/14.
# The scaling factor r_x is the ratio of the new width to the old width.
# r_x = (3/14) / 3 = 1/14.
r_x = 1/14

# The scaling factor in the y-direction, r_y, is (9/14) / 1 = 9/14.
# We don't need it for the final formula because the vertical fibers have dimension 0,
# as all black keys are at the same height.

# Step 2: Calculate the Minkowski-Bouligand dimension.
# The formula for this type of self-affine fractal is D = log(N) / log(1/r_x).
one_over_rx = 1 / r_x
dimension = math.log(N) / math.log(one_over_rx)

# Step 3: Print the results, including the formula and the numbers used.
print("The Minkowski-Bouligand dimension (D) for this fractal is calculated using the formula:")
print("D = log(N) / log(1/r_x)")
print("\nWhere:")
print(f"N (number of black keys) = {N}")
print(f"r_x (horizontal scaling factor) = {r_x:.4f}")
print("\nPlugging in the numbers:")
print(f"D = log({N}) / log(1/({r_x:.4f}))")
print(f"D = log({N}) / log({one_over_rx})")
print("\nResult:")
print(f"The dimension is: {dimension}")
