# The problem is to find the sum of the squares of the coefficients of the polynomial P(x).
# As derived in the plan, this sum is equal to the constant term of P(x)P(x^-1).
# The constant term of each factor in the product P(x)P(x^-1) is 4.
# There are 20 such factors (for i from 0 to 19).
# Therefore, the final result is 4 raised to the power of 20.

# Define the base and the exponent for the calculation.
base = 4
exponent = 20

# Calculate the result.
result = base ** exponent

# Print the final equation with the numbers involved.
print(f"The sum of the squares of the coefficients is given by the equation:")
print(f"{base} ** {exponent} = {result}")