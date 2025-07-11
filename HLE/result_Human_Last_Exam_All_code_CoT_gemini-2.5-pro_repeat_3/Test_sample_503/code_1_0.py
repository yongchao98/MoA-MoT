import math

# Step 1: Define the parameters from the problem.
# We are looking for adjunctions from [m] to [n].
m = 23
n = 37

# Step 2: The number of such adjunctions is given by the binomial coefficient
# C(m + n, n) or equivalently C(m + n, m).
# This is derived by counting the number of order-preserving maps R: [n] -> [m]
# that satisfy R(n) = m. This count is equivalent to the number of
# order-preserving maps from [n-1] to [m].
# Using the standard formula for counting order-preserving maps from [k] to [l],
# which is C(l + k + 1, k + 1), we set k = n-1 = 36 and l = m = 23.
# The number is C(23 + 36 + 1, 36 + 1) = C(60, 37).

# Step 3: Calculate the parameters for the binomial coefficient.
N = m + n
K = n

# Step 4: Calculate the binomial coefficient C(N, K).
result = math.comb(N, K)

# Step 5: Print the final equation and the result.
# The numbers in the final equation are N, K, and the result.
print(f"The number of internal adjunctions from [{m}] to [{n}] is calculated by the formula C({m} + {n}, {n}).")
print(f"C({N}, {K}) = {result}")