# A number is a "segmented number" if it cannot be written as the sum of
# two or more consecutive positive integers. This property is true for all
# powers of two (e.g., 1, 2, 4, 8, ...).

# The nth segmented number corresponds to the nth term in the sequence of powers of two,
# starting with 2**0. The formula for the nth term is 2**(n-1).

# We need to find the 50th element in this sequence.
n = 50

# The calculation is based on the formula 2**(n-1).
base = 2
exponent = n - 1

# Calculate the result
result = base ** exponent

# Print the final equation with all its numbers, as requested.
print(f"{base} ** {exponent} = {result}")