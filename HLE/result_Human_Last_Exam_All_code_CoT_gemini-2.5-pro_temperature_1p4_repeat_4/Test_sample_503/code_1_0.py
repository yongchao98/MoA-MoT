import math

# Define the parameters of the problem
# We are looking for adjunctions from [m] to [n]
m = 23
n = 37

# The number of such adjunctions is given by the binomial coefficient C(n+m, m).
# This is because an adjunction corresponds to an order-preserving map L:[m]->[n]
# with L(0)=0. Such a map is determined by a non-decreasing sequence
# L(1), ..., L(m) of m numbers chosen from {0, ..., n}.
# The number of such sequences is C((n+1)+m-1, m) = C(n+m, m).

# The total number in the 'n choose k' formula, C(N, K)
N = n + m
K = m

# Calculate the binomial coefficient C(N, K)
result = math.comb(N, K)

# Output the explanation and the final equation with all numbers
print(f"The number of internal adjunctions from [{m}] to [{n}] is given by the binomial coefficient C({n} + {m}, {m}).")
print(f"This is C({N}, {K}).")
print(f"C({N}, {K}) = {result}")
