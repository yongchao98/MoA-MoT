import math

# Step 1: Define the parameters from the problem description.

# N is the number of self-similar pieces. In a standard octave, there are 5 black keys.
# At each stage of the fractal construction, the keyboard structure is replicated on each of these 5 keys.
N = 5

# To find the scaling factors, we compare the dimensions of the whole keyboard
# with the dimensions of a single black key where the scaled copy is placed.

# The whole keyboard has a width of 3 units and a height of 1 unit.
keyboard_width = 3.0
keyboard_height = 1.0

# The problem states a black key is half as wide as a white key.
# With 7 equally sized white keys over a width of 3, one white key is 3/7 wide.
# So, a black key's width is (3/7) / 2 = 3/14.
black_key_width = 3.0 / 14.0

# A black key's height is 9/14 of the keyboard's height.
black_key_height = 9.0 / 14.0

# The scaling factors (r_x, r_y) are the ratios of the new dimensions to the old dimensions.
# r_x is the scaling factor in the x-direction (width).
r_x = black_key_width / keyboard_width

# r_y is the scaling factor in the y-direction (height).
r_y = black_key_height / keyboard_height


# Step 2: Explain and apply the box-counting dimension concept.
print("To find the Minkowski–Bouligand (box-counting) dimension, we analyze the self-affine structure of the object.")
print("The fractal is composed of N smaller copies of itself, scaled differently in width and height.\n")

print("The parameters of the fractal are:")
print(f"1. Number of self-similar copies (one for each black key): N = {N}")
print(f"2. Scaling factor in width: r_x = (3/14) / 3 = {r_x:.4f} (or 1/14)")
print(f"3. Scaling factor in height: r_y = (9/14) / 1 = {r_y:.4f} (or 9/14)\n")

# Step 3: Formulate the final equation.
# For self-affine fractals, the box-counting dimension D can be determined by analyzing
# how the number of boxes needed to cover the set grows as the box size shrinks.
# This analysis leads to a simplified formula D = log(A) / log(B).

# The base for the numerator, A, relates to the "effective" multiplication factor. It combines
# the number of new branches (N) with the ratio of the scaling factors, which accounts for the shape change.
numerator_base = N * (r_y / r_x)

# The base for the denominator, B, relates to the scale of the grid used for counting.
# We must use the smaller scaling factor (r_x) to ensure we capture the finest details of the fractal.
denominator_base = 1.0 / r_x

print("The box-counting dimension D is calculated using the formula: D = log(A) / log(B)")
print("where A and B are derived from the fractal's scaling properties.\n")

print("The base for the numerator of the formula, A, is:")
print(f"A = N * (r_y / r_x) = {N} * ( (9/14) / (1/14) ) = {N} * 9")
print(f"A = {int(numerator_base)}")

print("\nThe base for the denominator of the formula, B, is:")
print(f"B = 1 / r_x = 1 / (1/14)")
print(f"B = {int(denominator_base)}")

print("\nThus, the final equation for the dimension is:")
print(f"D = log({int(numerator_base)}) / log({int(denominator_base)})")

# Step 4: Calculate the final result.
dimension = math.log(numerator_base) / math.log(denominator_base)

print(f"\nThe calculated dimension is: {dimension}")