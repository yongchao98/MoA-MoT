# The problem asks for the 50th segmented number.
# Based on mathematical analysis, segmented numbers are powers of 2 (2^0, 2^1, 2^2, ...).
# The n-th segmented number is therefore 2^(n-1).

# Set the desired term number
n = 50

# The base of the power is 2
base = 2

# The exponent is n - 1
exponent = n - 1

# Calculate the result
result = base ** exponent

# The problem requires outputting each number in the final equation.
# The final equation is base^exponent = result.
print(f"{base}^{exponent} = {result}")