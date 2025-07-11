# The problem asks for the 50th segmented number.
# Based on mathematical analysis, segmented numbers are powers of 2 (1, 2, 4, 8, ...).
# The N-th segmented number is given by the formula 2^(N-1).

# For the 50th element, N = 50.
# The calculation is 2^(50-1) = 2^49.

# Define the components of the final equation
base = 2
n = 50
exponent = n - 1

# Calculate the result
result = base ** exponent

# Print the final equation with each number, as requested.
print(f"The equation to find the {n}th segmented number is: result = base ^ (n - 1)")
print(f"For the {n}th number:")
print(f"base = {base}")
print(f"n = {n}")
print(f"exponent = {n} - 1 = {exponent}")
print(f"result = {base} ^ {exponent} = {result}")