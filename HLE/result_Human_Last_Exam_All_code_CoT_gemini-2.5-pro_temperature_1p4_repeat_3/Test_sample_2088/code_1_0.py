# Based on the analytical solution, the value of the definite integral part of the expression
# is found to be (5**(1/4))/12.
# This script computes the final expression using this value.

# Define the numbers from the expression
base = 12
exponent = 4

# The value of the integral, I
# I = (5^(1/4)) / 12
integral_value = 5**(1/4) / 12

# The expression to compute is (base^exponent) * (I^exponent)
# This can be simplified to (base * I)^exponent
final_result = (base * integral_value)**exponent

# Output the components of the final calculation as an equation
print(f"The calculation is ({base})^{exponent} * ({integral_value})^{exponent}, which simplifies to ({base * integral_value})^{exponent}.")
print(f"Final result: {round(final_result)}")
