import math

# Based on the analytical derivation, the value of the integral expression is 5.
# Let's break down the final calculation as requested.
# The original expression is (12)^4 * (I_1 - I_2)^4
# We found that I = I_1 - I_2 = (5^(1/4)) / 12

# So, the expression becomes:
# (12)^4 * ( (5^(1/4)) / 12 )^4
# This simplifies to:
# 12^4 * ( 5 / 12^4 )
# Which is equal to 5.

base = 12
power = 4
numerator_of_integral_term = 5

# Calculate the result
result = (base ** power) * (numerator_of_integral_term / (base ** power))

# Print the final equation with each number.
# Here we show the simplified equation right before the final cancellation.
print("The expression simplifies to the following calculation:")
print(f"{base}^{power} * ({numerator_of_integral_term} / {base}^{power}) = {int(result)}")
