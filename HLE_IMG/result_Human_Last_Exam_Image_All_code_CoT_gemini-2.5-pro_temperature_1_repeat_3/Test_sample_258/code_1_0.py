import math

# The length of segment AB is given.
ab_length = 18468

# From the geometric derivation, we found that the length of CP is related to the length of AB
# by the golden ratio, φ (phi).
# The formula is: CP = AB * ( (sqrt(5) + 1) / 2 )

# Calculate the length of CP using the formula.
cp_length = ab_length * ((math.sqrt(5) + 1) / 2)

# As requested, here is the final equation with all the numbers.
print(f"The equation to find the length of CP is:")
print(f"CP = {ab_length} * ( (math.sqrt(5) + 1) / 2 )")

# And here is the final result.
print(f"The length of CP is {cp_length:.4f}")