# The problem asks for the 50th segmented number.
# A segmented number is a positive integer that cannot be expressed as the
# sum of two or more consecutive smaller positive integers.
# It can be shown mathematically that these numbers are the powers of 2:
# 1st: 2**0 = 1
# 2nd: 2**1 = 2
# 3rd: 2**2 = 4
# ...
# Therefore, the 50th segmented number is 2**(50-1), which is 2**49.

# The final equation is result = base ** exponent
n = 50
base = 2
exponent = n - 1
result = base ** exponent

# The final equation and its numbers (base, exponent, result) are printed below.
print(f"{base}**{exponent} = {result}")