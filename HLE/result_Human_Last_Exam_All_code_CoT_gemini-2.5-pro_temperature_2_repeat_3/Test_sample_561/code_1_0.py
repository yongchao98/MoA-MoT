import math

# Step 1: Define the parameters of the fractal system based on the problem description.
# N is the number of black keys in an octave.
N = 5

# s_x is the scaling factor in the horizontal (width) direction.
# Original width W = 3. Black key width w_b = 3/14.
# s_x = w_b / W = (3/14) / 3 = 1/14.
s_x = 1/14

# The scaling factor in the y-direction (s_y) is 9/14.
# We confirm the condition N * s_y > 1: 5 * (9/14) = 45/14 > 1.
# The formula for the dimension D is D = 1 + log(N) / log(1/s_x).

# Step 2: Calculate the terms for the formula.
# The inverse of the x-scaling factor.
inv_sx = 1 / s_x

# The numerator of the fractional part of the dimension formula.
log_N = math.log(N)

# The denominator of the fractional part of the dimension formula.
log_inv_sx = math.log(inv_sx)

# Step 3: Calculate the final dimension.
dimension = 1 + log_N / log_inv_sx

# Step 4: Print the full equation and the final result.
# The formula is D = 1 + log(N) / log(1/s_x)
# Note that log base is irrelevant as log(a)/log(b) is constant for any base.
print("The dimension D is calculated using the formula: D = 1 + log(N) / log(1/s_x)")
print("Plugging in the values:")
print(f"D = 1 + log({N}) / log({int(inv_sx)})")
print(f"D = {dimension}")
