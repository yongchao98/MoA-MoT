import math

# The problem is to find the number of internal adjunctions from [n] to [m].
# n and m are the parameters of the objects in the simplex category.
n = 23
m = 37

# The number of such adjunctions is given by the binomial coefficient C(n+m, n).
# This is derived from counting the number of order-preserving maps f: [n] -> [m]
# with the property f(0) = 0.
# For n=23 and m=37, we need to calculate C(23+37, 23) = C(60, 23).
N = n + m
K = n

# Calculate the binomial coefficient C(N, K)
result = math.comb(N, K)

# Print the final equation and the result.
# The user prompt requested to output the equation with each number.
# The equation is C(n+m, n) = result.
print(f"C({n} + {m}, {n}) = C({N}, {K}) = {result}")
