# The problem asks for the 50th segmented number.
# As established in the plan, the sequence of segmented numbers corresponds to powers of 2.
# The n-th segmented number is given by the formula 2^(n-1).

# Set the desired element number
n = 50

# The base for the power is 2
base = 2

# The exponent is n-1
exponent = n - 1

# Calculate the final result
result = base**exponent

# As requested, we print the final equation including each number.
# This shows that the 50th element is 2 raised to the power of (50-1), which is 49.
print(f"The 50th segmented number = {base}^({n}-1) = {result}")