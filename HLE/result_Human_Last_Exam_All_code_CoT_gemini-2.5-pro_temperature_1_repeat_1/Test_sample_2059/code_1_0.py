# The problem is to find the sum of the squares of the coefficients of the polynomial P(x).
# This sum is equal to the constant term of the product P(x) * P(1/x).
# As derived in the explanation, the constant term of each factor product is 4.
# There are 20 such factors in the product (for i from 0 to 19).
# Therefore, the final result is 4 raised to the power of 20.

# Base of the exponentiation
base = 4

# The exponent (number of terms in the product from i=0 to 19)
exponent = 20

# Calculate the result
result = base ** exponent

# Print the final equation with all numbers
print(f"The sum of the squares of the coefficients is calculated by the equation:")
print(f"{base}^{exponent} = {result}")
