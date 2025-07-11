# Based on the analytical derivation, the value of the integral is (5^(1/4))/12.
# We substitute this value into the expression.

# Define the numbers in the expression
base = 12
power = 4
integral_numerator_base = 5

# Calculate the terms of the equation
term1 = base ** power
term2_numerator = integral_numerator_base
term2_denominator = base ** power

# Calculate the final result
result = term1 * (term2_numerator / term2_denominator)

# Print the final equation with each number shown
print(f"The expression is ({base})^{power} * (integral)^4.")
print(f"The integral evaluates to ({integral_numerator_base}^(1/4)) / {base}.")
print(f"Plugging this into the expression gives:")
print(f"({base})^{power} * (({integral_numerator_base}^(1/4)) / {base})^{power} = ({term1}) * ({term2_numerator} / {term2_denominator}) = {int(result)}")