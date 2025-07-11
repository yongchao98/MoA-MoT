import math

# The final expression to be evaluated is (3/2) * 10^(10/3) + 37/4.
# We will calculate this value.
term1_coeff_num = 3
term1_coeff_den = 2
base = 10
exponent_num = 10
exponent_den = 3
term2_num = 37
term2_den = 4

# Calculate the result
result = (term1_coeff_num / term1_coeff_den) * (base**(exponent_num / exponent_den)) + (term2_num / term2_den)

# Print the final equation and its result
print(f"The final simplified equation is:")
print(f"({term1_coeff_num}/{term1_coeff_den}) * {base}^({exponent_num}/{exponent_den}) + ({term2_num}/{term2_den})")
print(f"\nThe calculated value is:")
print(result)