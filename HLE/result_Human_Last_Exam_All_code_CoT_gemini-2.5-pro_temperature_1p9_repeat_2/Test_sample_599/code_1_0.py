# The problem asks for the 50th segmented number.
# As explained in the plan, the sequence of segmented numbers corresponds to powers of 2.
# The n-th segmented number is given by the formula 2^(n-1).

# For this problem, we need the 50th element.
n = 50

# The base of the power is 2.
base = 2

# The exponent is n - 1.
exponent = n - 1

# Calculate the final result.
result = base**exponent

# Print the final equation showing the numbers involved.
print(f"The 50th segmented number is {base}^{exponent}, which is equal to:")
print(result)
