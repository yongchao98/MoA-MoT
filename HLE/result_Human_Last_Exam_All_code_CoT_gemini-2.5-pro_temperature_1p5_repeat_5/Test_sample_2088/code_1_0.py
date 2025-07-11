import math

# Define the constants from the expression
base = 12
exponent = 4

# From the analytical steps, we found the value of the definite integral part.
# The integral evaluates to (5^(1/4)) / 12.
integral_numerator_base = 5
integral_denominator = 12

# We are asked to compute (base^exponent) * (integral_value^exponent)
# Let's show the numbers in the final equation.

print("The integral part of the expression simplifies to:")
print(f"({integral_numerator_base}^(1/{exponent})) / {integral_denominator}")
print("\nNow, we compute the full expression:")
# This line shows the structure of the final calculation with all its numerical components.
print(f"({base})^{exponent} * (({integral_numerator_base}^(1/{exponent})) / {integral_denominator})^{exponent}")

# Step-by-step simplification of the final calculation
print(f"= ({base})^{exponent} * ({integral_numerator_base} / {integral_denominator}^{exponent})")
print(f"= {integral_numerator_base}")

# Compute the final result numerically for verification.
final_result = base**exponent * ((integral_numerator_base**(1/float(exponent))) / integral_denominator)**exponent
print(f"\nThe final answer is: {int(round(final_result))}")
