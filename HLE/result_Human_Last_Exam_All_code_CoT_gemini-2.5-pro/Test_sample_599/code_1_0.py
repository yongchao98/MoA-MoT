# The goal is to find the 50th segmented number.

# Step 1: Understand the definition of a segmented number.
# A positive integer N is a segmented number if it cannot be expressed as the
# sum of two or more consecutive smaller positive integers.

# Step 2: Derive the mathematical property of these numbers.
# A number N can be written as a sum of k >= 2 consecutive integers starting
# from a >= 1 if and only if N is not a power of 2.
# The reasoning is that such a sum can be written as N = k/2 * (2a + k - 1),
# which means 2N = k * (2a + k - 1). The factors k and (2a+k-1) have
# different parity, which implies 2N must have an odd factor greater than 1.
# Numbers whose double (2N) has no odd factor greater than 1 are powers of 2.
# Therefore, the numbers that *cannot* be written as such a sum are powers of 2.

# Step 3: Identify the sequence of segmented numbers.
# The segmented numbers are the powers of 2:
# 1st element: 2^0 = 1
# 2nd element: 2^1 = 2
# 3rd element: 2^2 = 4
# and so on.
# The n-th segmented number follows the formula 2^(n-1).

# Step 4: Calculate the 50th element using the formula.
# We need to find the 50th element, so we set n = 50.
# The 50th element is 2^(50 - 1) = 2^49.

# The target term in the sequence
n = 50

# The n-th segmented number is 2^(n-1)
exponent = n - 1
result = 2**exponent

# Print the final equation showing all the numbers involved in the calculation.
print(f"2 ** {exponent} = {result}")