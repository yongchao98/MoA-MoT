import math

# Step 5: Calculate the final result using the derived value of the integral.
# The expression is (12)^4 * I^4
# From our derivation, the integral I = (5^(1/4)) / 12.

coefficient = 12
power = 4
integral_numerator_base = 5

# The value of the integral I
# I = (integral_numerator_base ^ (1/4)) / coefficient
# We need to compute (coefficient^power) * (I^power)

# Let's compute I^4 first
I_pow_4_numerator = integral_numerator_base
I_pow_4_denominator = coefficient ** power

# Final expression
result = (coefficient ** power) * (I_pow_4_numerator / I_pow_4_denominator)

# Output the steps of the final calculation
print(f"The expression to compute is ({coefficient})^{power} * (I)^{power}")
print(f"We found the integral I = ({integral_numerator_base}^(1/4)) / {coefficient}")
print(f"So, we need to compute ({coefficient})^{power} * (({integral_numerator_base}^(1/4)) / {coefficient})^{power}")
print(f"= ({coefficient ** power}) * (({integral_numerator_base}) / ({coefficient ** power}))")
print(f"= {int(result)}")