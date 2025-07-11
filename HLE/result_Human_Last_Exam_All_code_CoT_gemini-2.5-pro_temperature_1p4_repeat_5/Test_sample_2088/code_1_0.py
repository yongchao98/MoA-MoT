import math

# Define the numbers from the problem
base = 12
power = 4
integral_result_numerator_base = 5

# The value of the integral is (5^(1/4)) / 12.
# We need to compute (12^4) * (integral_value^4).
# This is equivalent to (12 * integral_value)^4.

# Calculate the term inside the parentheses
term = base * (math.pow(integral_result_numerator_base, 1/power) / base)

# Raise the term to the power of 4
final_result = math.pow(term, power)

# Print the final equation with all the numbers
print(f"The calculation is: ({base})^({power}) * (({integral_result_numerator_base}^(1/{power})) / {base})^({power}) = ({base} * {integral_result_numerator_base}^(1/{power}) / {base})^({power}) = ({integral_result_numerator_base}^(1/{power}))^({power}) = {int(final_result)}")