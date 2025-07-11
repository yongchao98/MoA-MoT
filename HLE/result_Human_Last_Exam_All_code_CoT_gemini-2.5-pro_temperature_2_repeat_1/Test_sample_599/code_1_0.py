# The segmented numbers are powers of 2 (2^0, 2^1, 2^2, ...).
# The n-th segmented number is 2^(n-1).
# For the 50th element, n = 50.

# The equation is to find the value of 2^(50-1).
base = 2
n = 50
exponent = n - 1

# Calculate the result
result = base ** exponent

# Print the final equation and its result
# This line fulfills the requirement to "output each number in the final equation"
print(f"{base} ** {exponent} = {result}")