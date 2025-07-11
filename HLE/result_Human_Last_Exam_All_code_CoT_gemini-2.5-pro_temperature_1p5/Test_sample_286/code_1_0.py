import math

# Number of terms in the sequence
N = 100000

# The condition for M is M >= N * log10(2)
# We need to find the smallest integer M that satisfies this.

log10_2 = math.log10(2)
M_lower_bound = N * log10_2

# Since M must be an integer, we take the ceiling of the result.
M = math.ceil(M_lower_bound)

print(f"The number of terms is N = {N}")
print(f"The value of log10(2) is approximately: {log10_2}")
print(f"The lower bound for M is N * log10(2) = {N} * {log10_2} = {M_lower_bound}")
print(f"The smallest integer M must be the ceiling of this value.")
print(f"M = ceil({M_lower_bound}) = {M}")
