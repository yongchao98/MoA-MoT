# The problem asks for the 50th segmented number.
# A segmented number is a positive integer that cannot be expressed as the sum of
# two or more consecutive positive integers.
# This property holds true for all powers of 2.
# The sequence of segmented numbers is: 2^0, 2^1, 2^2, 2^3, ...
# So, the nth segmented number is 2^(n-1).

# We need to find the 50th element, so n = 50.
n = 50

# The formula is 2^(n-1).
base = 2
exponent = n - 1

# Calculate the result.
result = base ** exponent

# Print the final equation with all its numbers.
# The numbers in the equation are the base (2), the exponent (49), and the result.
print(f"The {n}th segmented number is given by the equation:")
print(f"{base} ** {exponent} = {result}")
