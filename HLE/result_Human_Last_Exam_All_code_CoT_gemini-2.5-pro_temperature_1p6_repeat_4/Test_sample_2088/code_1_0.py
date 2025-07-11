# The original expression is (12)^4 * (I)^4
# From the derivation, the value of the integral I is (5^(1/4)) / 12.

# Constants from the expression
c1 = 12
integral_val_numerator_base = 5
integral_val_denominator = 12

# Calculate the first term: (12)^4
term1 = c1**4

# The integral value is integral_val = (5^(1/4)) / 12
# We need to compute (integral_val)^4
integral_val_pow_4_numerator = integral_val_numerator_base**(0.25 * 4)
integral_val_pow_4_denominator = integral_val_denominator**4

# The second term of the expression
term2 = integral_val_pow_4_numerator / integral_val_pow_4_denominator

# The final result is the product of the two terms
result = term1 * term2

# Print the calculation steps
print(f"The expression is ({c1})^4 * (({integral_val_numerator_base}^(1/4)) / {integral_val_denominator})^4")
print(f"= {term1} * (({integral_val_numerator_base}^1) / {integral_val_denominator**4})")
print(f"= {term1} * ({int(integral_val_pow_4_numerator)} / {integral_val_pow_4_denominator})")
print(f"= {int(result)}")