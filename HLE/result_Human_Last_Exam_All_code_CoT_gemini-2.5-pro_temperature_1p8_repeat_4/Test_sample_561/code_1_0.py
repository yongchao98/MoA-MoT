import math

# Step 1: Define the parameters based on the problem description.
# N is the number of black keys in an octave.
N = 5

# The keyboard's overall dimensions.
keyboard_width = 3.0
keyboard_height = 1.0

# The properties of the black keys.
black_key_height_ratio = 9.0 / 14.0
# A black key is half the width of a white key. There are 7 white keys.
black_key_width = (keyboard_width / 7.0) / 2.0

# Step 2: Calculate the horizontal and vertical scaling factors.
# s_w is the ratio of a black key's width to the total keyboard's width.
s_w = black_key_width / keyboard_width
# s_h is the ratio of a black key's height to the total keyboard's height.
s_h = black_key_height_ratio * keyboard_height / keyboard_height

# Step 3: Define the terms for the dimension formula.
# The formula for the dimension D is: log(N * s_h / s_w) / log(1 / s_w)
numerator_arg = N * (s_h / s_w)
denominator_arg = 1 / s_w

# Step 4: Calculate the dimension.
dimension = math.log(numerator_arg) / math.log(denominator_arg)

# Step 5: Print the explanation and the final result.
print("To find the Minkowski–Bouligand dimension of the black keys, we model it as a self-affine fractal.")
print(f"The number of self-similar pieces is N = {N} (the number of black keys in an octave).")
print(f"The horizontal scaling factor is s_w = {s_w:.5f}")
print(f"The vertical scaling factor is s_h = {s_h:.5f}\n")

print("The dimension D is found using the formula:")
print("D = log(N * s_h / s_w) / log(1 / s_w)\n")

print("Plugging in the numbers, the final equation is:")
print(f"D = log({int(numerator_arg)}) / log({int(denominator_arg)})")

print(f"\nThe calculated Minkowski–Bouligand dimension is: {dimension}")
