# The problem is to find the 50th segmented number.
# Through mathematical analysis, it's determined that segmented numbers
# are powers of 2 (i.e., numbers of the form 2^k).
# The first segmented number is 1 (2^0).
# The second is 2 (2^1).
# The nth segmented number is 2^(n-1).

# We need to find the 50th element, so n = 50.
n = 50

# The corresponding exponent is n - 1.
exponent = n - 1
base = 2

# Calculate the result of 2^49.
result = base ** exponent

# Per the instruction to "output each number in the final equation",
# the following line prints the base (2), the exponent (49), and the result.
print(f"{base} ** {exponent} = {result}")