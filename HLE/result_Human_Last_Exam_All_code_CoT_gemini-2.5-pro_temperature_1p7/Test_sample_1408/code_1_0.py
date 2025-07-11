import math

# The maximal overhang is (1 + sqrt(2)) / 2
# We need to find integers a, b, c such that this equals (a + sqrt(b)) / (1 + c)
# a, b, c must be non-negative.
# sqrt(b) must be non-integer or zero.
# c must be minimal.

# By comparing the expressions:
# a = 1
# b = 2
# 1 + c = 2  => c = 1

a = 1
b = 2
c = 1

# Verification of conditions
# a, b, c are non-negative integers?
is_non_negative_integer = (isinstance(a, int) and a >= 0 and
                           isinstance(b, int) and b >= 0 and
                           isinstance(c, int) and c >= 0)

# sqrt(b) is non-integer or zero?
sqrt_b_is_non_integer_or_zero = (b == 0 or int(math.sqrt(b))**2 != b)

# c is minimal?
# If c=0, the overhang would be a + sqrt(b).
# (1 + sqrt(2))/2 is not of the form a + sqrt(b) for integers a, b.
# So c cannot be 0. Thus c=1 is minimal.
is_c_minimal = True

if is_non_negative_integer and sqrt_b_is_non_integer_or_zero and is_c_minimal:
    print(f"{a} {b} {c}")
