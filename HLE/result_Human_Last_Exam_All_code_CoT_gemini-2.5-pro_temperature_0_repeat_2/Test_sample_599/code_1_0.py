# The goal is to find the 50th segmented number.
# Step 1: Understand that segmented numbers are powers of 2.
# The sequence is 2^0, 2^1, 2^2, ... which is 1, 2, 4, ...
# The n-th segmented number is 2^(n-1).

# Step 2: Calculate the 50th element.
# For n=50, the segmented number is 2^(50-1) = 2^49.

# The position in the sequence.
n = 50

# The base of the power.
base = 2

# The exponent for the n-th element.
exponent = n - 1

# Calculate the result.
result = base ** exponent

# Step 3: Output the numbers in the final equation as requested.
# The final equation is base^exponent = result.
print(f"The 50th segmented number is found from the equation: {base} ** {exponent}")
print(f"The result is: {result}")
print(f"So, the final equation is: {base} ** {exponent} = {result}")
