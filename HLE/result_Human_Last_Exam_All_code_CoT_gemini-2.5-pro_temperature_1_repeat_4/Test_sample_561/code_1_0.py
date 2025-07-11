import math

# Step 1: Define the parameters of the fractal construction.
# N is the number of self-similar copies. In a piano octave, there are 5 black keys.
N = 5

# Step 2: Calculate the scaling factors.
# The original keyboard is 3 units wide and 1 unit high.
keyboard_width = 3.0
keyboard_height = 1.0

# The black keys are 9/14 of the height.
black_key_height = (9.0 / 14.0) * keyboard_height

# The 7 white keys make up the width, and black keys are half as wide as white keys.
white_key_width = keyboard_width / 7.0
black_key_width = white_key_width / 2.0

# The scaling factors are the ratio of the black key dimensions to the keyboard dimensions.
r_x = black_key_width / keyboard_width
r_y = black_key_height / keyboard_height

# The singular values are the scaling factors, sorted from largest to smallest.
sigma_1 = max(r_x, r_y)
sigma_2 = min(r_x, r_y)

# Step 3: Define the equation for the dimension D.
# For a self-affine set where 1 <= D < 2, the dimension D is the solution to:
# N * sigma_1 * sigma_2^(D-1) = 1
# Substituting the numbers: 5 * (9/14) * (1/14)^(D-1) = 1
# This simplifies to the form: 45/14 * (1/14)^(D-1) = 1
# Taking logs and rearranging gives: D = log(45) / log(14)

# We define the numerator and denominator for the final simplified equation D = log(A)/log(B)
# A = N * sigma_1 / sigma_2 = 5 * (9/14) / (1/14) = 5 * 9 = 45
# B = 1 / sigma_2 = 1 / (1/14) = 14
A = 45
B = 14

# Step 4: Calculate the dimension and print the results.
dimension = math.log(A) / math.log(B)

print("The Minkowski–Bouligand dimension (D) of this fractal is calculated for a self-affine set.")
print("The dimension D is the solution to the equation: N * \u03C3\u2081 * \u03C3\u2082^(D-1) = 1")
print(f"Number of copies (N): {N}")
print(f"Largest scaling factor (\u03C3\u2081): {sigma_1:.4f} (which is 9/14)")
print(f"Smallest scaling factor (\u03C3\u2082): {sigma_2:.4f} (which is 1/14)")
print("\nSubstituting the values, the equation is:")
print(f"{N} * (9/14) * (1/14)^(D-1) = 1")
print("\nThis equation simplifies to D = log(A) / log(B)")
print(f"Where A = {A} and B = {B}")
print(f"\nFinal Equation: D = log({A}) / log({B})")
print(f"\nThe calculated dimension is: {dimension}")