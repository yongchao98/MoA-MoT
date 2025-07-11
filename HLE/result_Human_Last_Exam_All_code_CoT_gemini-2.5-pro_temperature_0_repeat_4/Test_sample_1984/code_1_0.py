# The number of stable equilibrium points.
m = 1

# The exponent in the expression.
exponent = 4048

# Calculate 2 to the power of the exponent.
# Python handles arbitrarily large integers automatically.
power_of_two = 2**exponent

# Calculate the final result.
result = m - power_of_two

# Print the final equation and its result.
# The result is a very large negative number.
print(f"{m} - 2^{exponent} = {result}")