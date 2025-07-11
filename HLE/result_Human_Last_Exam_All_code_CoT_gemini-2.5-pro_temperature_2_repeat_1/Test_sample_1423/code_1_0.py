# The number of distinct digits allowed.
# To find the maximum possible number of digits in N, we should use the
# maximum number of distinct digits available.
m = 5

# The formula for the maximum length f(m) of a valid sequence with m distinct digits
# is f(m) = 2^m - 1.
base = 2

# Perform the calculation
power_val = base**m
result = power_val - 1

# Print the explanation and the result as a step-by-step equation.
print(f"Let m be the number of distinct digits. The problem states m <= 5.")
print(f"The maximum length of the sequence is given by the formula: 2^m - 1.")
print(f"To maximize the length, we choose m = {m}.")
print(f"The calculation is as follows:")
print(f"{base}^{m} - 1 = {power_val} - 1 = {result}")
