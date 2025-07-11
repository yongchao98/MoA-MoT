# 1. The problem asks for the 50th segmented number.
# 2. A segmented number is a positive integer that cannot be expressed as the sum of two or more consecutive positive integers.
# 3. It can be proven that the set of such numbers is the set of powers of 2: {2^0, 2^1, 2^2, ...}.
# 4. The first element is 1 (2^0), the second is 2 (2^1), and so on.
# 5. The n-th segmented number is given by the formula 2^(n-1).
# 6. To find the 50th element, we need to calculate 2^(50-1).

# The position in the sequence
n = 50

# The exponent is n-1
exponent = n - 1

# Calculate the result, which is 2 to the power of the exponent
result = 2**exponent

# Print the final equation showing all the numbers involved in the calculation
print(f"The {n}th segmented number is calculated by the equation: 2^({n} - 1) = 2^{exponent} = {result}")