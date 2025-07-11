import math

# Step 1: Define the parameters of the fractal based on the problem description.

# N is the number of self-affine pieces, which is the number of black keys in one octave.
N = 5

# s_x and s_y are the scaling factors in the x and y directions.
# A miniature keyboard (3x1) is scaled to fit into a black key (3/14 x 9/14).
# s_x = (width of black key) / (width of keyboard) = (3/14) / 3 = 1/14
# s_y = (height of black key) / (height of keyboard) = (9/14) / 1 = 9/14
s_x = 1/14
s_y = 9/14

# Step 2: Identify the larger and smaller scaling factors.
s_large = s_y
s_small = s_x

# Step 3: Use the formula for the box-counting dimension of this self-affine set:
# D = (log(N) + log(s_large / s_small)) / log(1 / s_small)
# We can simplify this to D = log(N * (s_large / s_small)) / log(1 / s_small)

# Calculate the arguments for the log functions in the simplified formula.
s_ratio = s_large / s_small
numerator_arg = N * s_ratio
denominator_arg = 1 / s_small

# Step 4: Calculate the final dimension.
dimension = math.log(numerator_arg) / math.log(denominator_arg)

# Step 5: Print the breakdown of the final equation and the result.
print("The Minkowski–Bouligand dimension D is calculated using the formula for self-affine fractals.")
print("The final equation for the dimension is:")
print("")
print(f"D = log({int(N)} * ({s_y:.2f} / {s_x:.2f})) / log(1 / {s_x:.2f})")
print(f"D = log({int(N)} * {int(s_ratio)}) / log({int(denominator_arg)})")
print(f"D = log({int(numerator_arg)}) / log({int(denominator_arg)})")
print("")
print("The calculated dimension is:")
print(dimension)
