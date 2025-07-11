import math
from fractions import Fraction

# Step 1: Define parameters from the problem description
# N is the number of self-similar copies, which equals the number of black keys in one octave.
N = 5
keyboard_width = 3
keyboard_height = 1

# Calculate the dimensions of a single black key based on the problem's description.
# A one-octave keyboard has 7 white keys of equal width.
white_key_width = Fraction(keyboard_width, 7)
# Black keys are half as wide as white keys.
black_key_width = white_key_width / 2
# Black keys cover 9/14 of the keyboard height.
black_key_height = Fraction(9, 14) * keyboard_height

# Step 2: Calculate the anisotropic scaling factors
# s_x is the ratio of the new width (black_key_width) to the old width (keyboard_width).
s_x = black_key_width / keyboard_width
# s_y is the ratio of the new height (black_key_height) to the old height (keyboard_height).
s_y = black_key_height / keyboard_height

print(f"The fractal is constructed from N = {N} copies per iteration.")
print(f"The scaling factor for the width is s_x = {s_x.numerator}/{s_x.denominator}.")
print(f"The scaling factor for the height is s_y = {s_y.numerator}/{s_y.denominator}.")
print("-" * 20)

# Step 3: Set up parameters for the dimension formula.
# The scaling factors must be sorted in descending order.
s_factors = sorted([s_x, s_y], reverse=True)
s1 = s_factors[0]
s2 = s_factors[1]
# For these parameters, the integer 'p' in the dimension formula is 1.
p = 1

print("Using the box-counting dimension formula for self-affine sets:")
print("D = p + log(N * s1) / log(1/s2)")
print(f"With N={N}, s1={s1}, s2={s2}, and p={p}, the formula is initially:")
print(f"D = 1 + log({N} * {s1}) / log(1 / {s2})")


# Step 4: Simplify the formula to its final elegant form.
# The formula D = 1 + log(45/14) / log(14) simplifies to D = log(45)/log(14).
# The numbers in the simplified formula are 45 and 14.
final_num = N * s1.numerator
final_den = int(1/s2)

# Step 5: Calculate the numerical result from the simplified formula.
dimension = math.log(final_num) / math.log(final_den)

print("\nAfter algebraic simplification, the final equation for the dimension (D) is:")
print(f"D = log({final_num}) / log({final_den})")
print(f"\nThe calculated dimension is: {dimension}")
