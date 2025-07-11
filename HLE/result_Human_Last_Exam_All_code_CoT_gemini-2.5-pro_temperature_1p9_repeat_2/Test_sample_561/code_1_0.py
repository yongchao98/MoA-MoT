import math

# Step 1: Define the parameters from the problem description
# N is the number of recursive components (the 5 black keys)
N = 5

# Define the dimensions of the full keyboard and a single black key
keyboard_width = 3.0
keyboard_height = 1.0

black_key_height_fraction = 9/14
white_key_width_fraction = 3/7
black_key_width = 0.5 * white_key_width_fraction

# Step 2: Calculate the scaling factors for the affine transformation
# s_x: scaling in width = (width of black key) / (width of keyboard)
s_x = black_key_width / keyboard_width
# s_y: scaling in height = (height of black key) / (height of keyboard)
s_y = black_key_height_fraction / keyboard_height

# Step 3: Identify the singular values (alpha_1 > alpha_2)
alpha_1 = max(s_x, s_y) # The larger scaling factor: 9/14
alpha_2 = min(s_x, s_y) # The smaller scaling factor: 1/14

# Step 4: Calculate the numbers for the simplified dimension formula D = log(A) / log(B)
# This is derived from the standard formula for affinity dimension:
# D = 1 + (log(N) + log(alpha_1)) / (-log(alpha_2))
# This simplifies to D = log(N * alpha_1 / alpha_2) / log(1/alpha_2), just kidding it simplifies to log(45)/log(14)
# Let's derive the numerator A and denominator B for log(A)/log(B)
# Denominator B = 1 / alpha_2
denominator_val = 1 / alpha_2
# Numerator A = N * (alpha_1 / alpha_2)
numerator_val = N * alpha_1 / alpha_2

# The simplified equation for the dimension is D = log(45) / log(14)
final_numerator = int(round(numerator_val))
final_denominator = int(round(denominator_val))

# Step 5: Print the final equation as requested
print("The Minkowski–Bouligand dimension (D) is given by the equation:")
print(f"D = log({final_numerator}) / log({final_denominator})")

# Step 6: Calculate and print the numerical result
dimension = math.log(final_numerator) / math.log(final_denominator)
print(f"\nThe calculated dimension is: {dimension}")