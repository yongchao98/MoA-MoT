import math
from fractions import Fraction

# Step 1: Define the parameters from the problem description.
# Number of self-similar copies (number of black keys in an octave).
N = 5

# Dimensions of the main keyboard.
keyboard_w = 3
keyboard_h = 1

# Dimensions of the black keys.
# There are 7 white keys making up the width.
white_key_w = Fraction(keyboard_w, 7)
# Black keys are half as wide as white keys.
black_key_w = white_key_w / 2
# Black keys cover 9/14 of the height.
black_key_h = Fraction(9, 14) * keyboard_h

# Step 2: Calculate the scaling factors rx and ry.
# rx is the ratio of a black key's width to the keyboard's width.
rx = black_key_w / keyboard_w
# ry is the ratio of a black key's height to the keyboard's height.
ry = black_key_h / keyboard_h

# Step 3: Identify the singular values s1 and s2 for the dimension formula.
# s1 is the larger scaling factor, s2 is the smaller one.
s1 = max(rx, ry)
s2 = min(rx, ry)

# Step 4: Define the numbers for the final equation D = log(45) / log(14)
# The number 45 comes from N * s1 * (denominator of s2).
# 5 * (9/14) * 14 = 45
num_for_log = N * s1.numerator * s2.denominator / s1.denominator
# The number 14 is the denominator of the scaling factors.
den_for_log = s2.denominator

# Step 5: Calculate the dimension D.
D = math.log(num_for_log) / math.log(den_for_log)

# Step 6: Print the explanation and the result.
print("Calculating the Minkowski–Bouligand Dimension (D) for the fractal piano keys.")
print("-" * 60)
print(f"1. Number of self-similar copies (N): {N}")
print(f"2. Scaling factor in x-direction (rx = black_key_width / keyboard_width): {rx}")
print(f"3. Scaling factor in y-direction (ry = black_key_height / keyboard_height): {ry}")
print(f"4. The larger scaling factor (s1) is {s1}, and the smaller (s2) is {s2}.")
print("\nFor a self-affine fractal, the dimension D is the solution to:")
print("N * s1 * s2^(D-1) = 1")
print(f"Plugging in the numbers: {N} * ({s1}) * ({s2})^(D-1) = 1")
print("\nSolving this equation for D gives the simplified formula:")
print(f"D = log({int(num_for_log)}) / log({den_for_log})")
print("-" * 60)
print(f"The calculated dimension is D = {D}")
