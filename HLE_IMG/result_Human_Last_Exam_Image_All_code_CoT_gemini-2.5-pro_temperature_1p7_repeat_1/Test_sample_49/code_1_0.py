import math

# Step 1: Define the components symbolically.
# We are asked to express the cutoff frequency in terms of 'r' and 'C'.

# Step 2: The Thevenin resistance R_th, as derived from analyzing the
# infinite ladder network, is r * (1 + sqrt(3)).
# We calculate the numerical coefficient.
sqrt_3 = math.sqrt(3)
coefficient = 1 + sqrt_3

# Step 3: The cutoff angular frequency is w_c = 1 / (R_th * C).
# We will print the final equation showing the breakdown of the numbers involved,
# as requested by the user prompt.
one = 1.0

print("The cutoff angular frequency w_c is determined by the equation:")
print("w_c = 1 / (R_th * C)")
print("where the Thevenin resistance R_th = (1 + sqrt(3)) * r")
print("\nSubstituting the numerical value of sqrt(3):")
print(f"w_c = {one} / (({one} + {sqrt_3:.5f}) * r * C)")
print("\nThis simplifies to the final expression:")
print(f"w_c = {one} / ({coefficient:.5f} * r * C)")