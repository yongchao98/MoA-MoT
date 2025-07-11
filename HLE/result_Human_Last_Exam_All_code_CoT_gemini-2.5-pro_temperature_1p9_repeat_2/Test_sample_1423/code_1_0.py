# The number of distinct digits we can use. To maximize the length,
# we should use the maximum number allowed.
k = 5

# The base for the exponential growth in the formula L(k) = 2^k - 1
base = 2

# The formula for the maximum length of a sequence with k distinct digits
# that satisfies the given condition is 2^k - 1.
# Here we calculate the result for k=5.
result = base**k - 1

# As requested, output the numbers in the final equation.
# This shows the equation 2**5 - 1 = 31.
power_val = base**k
subtrahend = 1
print(f"{base}**{k} - {subtrahend} = {power_val} - {subtrahend} = {result}")