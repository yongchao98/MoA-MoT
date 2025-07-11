import math

# The problem is to count the number of internal adjunctions from [23] to [37].
# This is equivalent to counting the number of order-preserving maps f: [23] -> [37]
# such that f(0) = 0.
# This corresponds to choosing 23 numbers (f(1), ..., f(23)) from the range [0, 37]
# with replacement, in non-decreasing order.
# The number of ways to do this is given by the binomial coefficient C(n+k-1, k),
# where n is the number of items to choose from (38, i.e., 0 to 37) and k is the number
# of items to choose (23, i.e., f(1) to f(23)).
# So, we need to calculate C(38+23-1, 23) = C(60, 23).

n = 60
k = 23

# Calculate the binomial coefficient C(n, k)
result = math.comb(n, k)

# Print the final equation and the result
print(f"The number of internal adjunctions is given by the binomial coefficient C({n}, {k}).")
print(f"C({n}, {k}) = {result}")
