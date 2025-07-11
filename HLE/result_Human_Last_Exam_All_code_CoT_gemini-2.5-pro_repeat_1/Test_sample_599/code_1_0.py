# The problem asks for the 50th segmented number.
# A segmented number is a positive integer that cannot be expressed as the sum of
# two or more consecutive smaller positive integers.
# Through mathematical analysis, it can be shown that these numbers are the powers of 2.
# 1st segmented number = 2^(1-1) = 2^0 = 1
# 2nd segmented number = 2^(2-1) = 2^1 = 2
# ...
# Nth segmented number = 2^(N-1)

# We need to compute the 50th element of this sequence.
n = 50

# The formula for the nth element is 2^(n-1)
base = 2
exponent = n - 1

# Calculate the result
result = base ** exponent

# Print the final equation with all its components, as requested.
# The equation is: base ^ exponent = result
print(f"{base} ** {exponent} = {result}")