# The problem is to find the maximum possible number of digits in an integer N
# that uses at most 5 distinct digits, subject to a specific substring property.
# Based on the analysis, the maximum length L for a number with k distinct digits
# follows the formula L(k) = 2^k - 1.

# We want to find the maximum length for k <= 5. Since L(k) increases with k,
# the maximum will be at k=5.
k = 5

# The final equation to calculate the length is 2^k - 1.
base = 2
exponent = k
subtrahend = 1

# Calculate the result
result = base**exponent - subtrahend

# The final instruction asks to output each number in the final equation.
# The final equation is: 2^5 - 1 = 31
# The numbers are 2, 5, 1, and 31.
print(f"The calculation is based on the formula: {base}^{exponent} - {subtrahend}")
print(f"The result is: {result}")
print("The numbers that make up the final equation are:")
print(base)
print(exponent)
print(subtrahend)
print(result)
