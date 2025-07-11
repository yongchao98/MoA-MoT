# A number is a "segmented number" if it is a positive integer that cannot be
# written as the sum of two or more consecutive smaller positive integers.
# It's a known mathematical property that such numbers are exactly the powers of 2.
#
# The sequence of segmented numbers is:
# 1st element = 1 = 2^0
# 2nd element = 2 = 2^1
# 3rd element = 4 = 2^2
# nth element = 2^(n-1)
#
# We need to find the 50th element of this sequence.

# The position in the sequence
n = 50

# The base for our calculation is 2
base = 2

# The exponent is n - 1
exponent = n - 1

# Calculate the result
result = base ** exponent

# Print the final equation with all its components
print(f"The 50th segmented number is calculated as {base} raised to the power of {exponent}.")
print(f"{base} ** {exponent} = {result}")