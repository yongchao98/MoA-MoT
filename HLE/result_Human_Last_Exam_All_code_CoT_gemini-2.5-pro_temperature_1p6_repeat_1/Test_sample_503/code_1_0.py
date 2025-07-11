import math

# We need to calculate the number of internal adjunctions from [23] to [37].
# This is equivalent to the number of order-preserving maps g: [37] -> [23]
# such that g(37) = 23.
# This in turn is equivalent to the number of order-preserving maps from [36] to [23].
# The formula for the number of order-preserving maps from [n] to [m] is C(m+n+1, n+1).
# Here, n = 36 and m = 23.
n = 36
m = 23
N = m + n + 1
K = n + 1

# The number is C(60, 37)
result = math.comb(N, K)

print(f"The number of adjunctions is given by the binomial coefficient C({N}, {K}).")
print(f"C({N}, {K}) = {result}")