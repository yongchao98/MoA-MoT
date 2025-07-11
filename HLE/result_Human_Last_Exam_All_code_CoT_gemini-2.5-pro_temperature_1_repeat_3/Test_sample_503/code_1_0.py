import math

# The problem is to count the number of order-preserving maps L: [23] -> [37]
# with the property L(0) = 0.
# This is equivalent to counting the number of non-decreasing sequences
# of length 23, with elements from {0, 1, ..., 37}.

# Let k be the length of the sequence and N be the maximum value.
# k = 23
# N = 37
# The number of such sequences is given by the binomial coefficient C(N+k, k).

k = 23
N = 37

n_comb = N + k
k_comb = k

# Calculate the binomial coefficient C(60, 23)
result = math.comb(n_comb, k_comb)

# Print the final equation and the result
print(f"The number of internal adjunctions is C({n_comb}, {k_comb}) = {result}")