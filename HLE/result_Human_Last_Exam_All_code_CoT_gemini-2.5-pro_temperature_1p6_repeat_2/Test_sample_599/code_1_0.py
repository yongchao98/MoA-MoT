# The problem asks for the 50th segmented number.
# Based on the definition, a segmented number is a positive integer that cannot be
# expressed as the sum of two or more consecutive positive integers.
# A key result from number theory states that these numbers are precisely the powers of 2.
# The sequence of segmented numbers is:
# 1st term: 2^0 = 1
# 2nd term: 2^1 = 2
# 3rd term: 2^2 = 4
# and so on.
# The n-th term in the sequence is given by the formula 2^(n-1).

# We need to find the 50th term, so we set n = 50.
n = 50

# The base for our calculation is 2.
base = 2

# The exponent is n - 1.
exponent = n - 1

# Calculate the final result.
result = base ** exponent

# Print the result to the console.
print(result)