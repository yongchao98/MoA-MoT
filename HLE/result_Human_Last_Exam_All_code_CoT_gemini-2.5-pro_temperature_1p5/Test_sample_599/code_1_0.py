# The segmented numbers are powers of 2 (2^0, 2^1, 2^2, ...).
# The n-th segmented number is given by the formula 2^(n-1).
# We need to find the 50th segmented number.

# The position in the sequence.
n = 50

# The base for our calculation is 2.
base = 2

# The exponent is n - 1.
power = n - 1

# Calculate the result. Python's integers handle arbitrary size, so this is safe.
result = base**power

# As requested, print the final equation including each number.
print(f"The 50th segmented number is calculated by the equation:")
print(f"{base} ** {power} = {result}")