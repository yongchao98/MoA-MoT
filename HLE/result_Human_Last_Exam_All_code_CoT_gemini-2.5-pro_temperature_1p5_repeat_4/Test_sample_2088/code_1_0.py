# The original expression is (12)^4 * I^4
# where I is the value of the definite integral.

# From the analytical solution, we found the value of the integral I.
# I = (5^(1/4)) / 12

# Let's define the components of the original expression.
base_of_first_term = 12
power_of_first_term = 4
integral_value_numerator = 5
integral_value_denominator = 12
power_of_integral = 4

# Calculate the first term
term1 = base_of_first_term ** power_of_first_term

# Calculate the value of the integral raised to its power
powered_integral_numerator = integral_value_numerator
powered_integral_denominator = integral_value_denominator ** power_of_integral

# Calculate the final result
final_result = term1 * (powered_integral_numerator / powered_integral_denominator)

# Print out the equation with the calculated numbers
print(f"The calculation is: ({base_of_first_term})^{power_of_first_term} * ({integral_value_numerator} / {powered_integral_denominator})")
print(f"= {term1} * ({powered_integral_numerator} / {powered_integral_denominator})")
print(f"= {int(final_result)}")

# Also printing the equation in the symbolic form before final calculation
print("\nSymbolic calculation:")
print(f"({base_of_first_term})^{power_of_first_term} * ( ( {integral_value_numerator}^(1/4) / {integral_value_denominator} )^{power_of_integral} ) = {int(final_result)}")