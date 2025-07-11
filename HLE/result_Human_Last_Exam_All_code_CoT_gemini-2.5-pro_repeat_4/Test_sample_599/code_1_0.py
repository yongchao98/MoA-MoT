# The problem asks for the 50th segmented number.
# As derived in the plan, segmented numbers are powers of 2 (2^0, 2^1, 2^2, ...).
# The n-th segmented number follows the formula 2^(n-1).

# We need to find the 50th element.
n = 50

# The base is 2 for powers of two.
base = 2

# The exponent is n - 1.
exponent = n - 1

# Calculate the result. Python handles large integers automatically.
result = base ** exponent

# The prompt requires printing each number in the final equation.
# The following line prints the base, exponent, and the final result.
print(f"{base} ** {exponent} = {result}")