# The problem asks for the 50th segmented number.
# A segmented number is a positive integer that cannot be expressed as the sum of two or more consecutive smaller positive integers.
# It's a mathematical property that such numbers are precisely the powers of 2.
# 1st segmented number: 1 = 2^0 = 2^(1-1)
# 2nd segmented number: 2 = 2^1 = 2^(2-1)
# 3rd segmented number: 4 = 2^2 = 2^(3-1)
# Following this pattern, the n-th segmented number is 2^(n-1).
# Therefore, the 50th segmented number is 2^(50-1).

# Define the base and exponent for the calculation.
base = 2
exponent = 50 - 1

# Calculate the result.
result = base ** exponent

# Print the final equation with all its numbers.
print(f"The 50th segmented number is given by the equation:")
print(f"{base} ** {exponent} = {result}")
