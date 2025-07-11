# The segmented numbers are powers of 2 (2^0, 2^1, 2^2, ...).
# To find the 50th element in this sequence, we calculate 2 to the power of (50-1).

# The number of the term in the sequence
n = 50

# The base of the power is 2
base = 2

# The exponent is n - 1
exponent = n - 1

# Calculate the result
result = base ** exponent

# The final equation includes the base, exponent, and the result.
print(f"{base} ^ {exponent} = {result}")
