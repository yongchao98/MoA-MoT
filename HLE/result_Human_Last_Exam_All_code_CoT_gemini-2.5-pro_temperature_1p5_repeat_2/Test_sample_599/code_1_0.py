# The problem is to find the 50th "segmented number".
# A segmented number is a positive integer that cannot be written as the sum
# of two or more consecutive smaller positive integers.
#
# Our analysis shows that these numbers are the powers of 2.
# The sequence starts:
# 1st term: 1 = 2^0
# 2nd term: 2 = 2^1
# 3rd term: 4 = 2^2
# ...
# n-th term: 2^(n-1)
#
# We need to find the 50th term in this sequence.

# The term number we are looking for.
n = 50

# The base of the power calculation for segmented numbers.
base = 2

# The exponent for the n-th term is n-1.
exponent = n - 1

# Calculate the 50th segmented number.
# Python's integers have arbitrary precision, so they can handle this large number.
result = base ** exponent

# The final calculation is result = base ** exponent.
# As requested, we will output each number in this final equation.
print(f"The 50th segmented number is the result of {base} raised to the power of {exponent}.")
print(f"Final Equation: {base} ** {exponent} = {result}")
