import math

# Number of steps in the sequence
N = 100000

# The condition is L_max <= 2*S, where S = 1 - 10**(-M)
# L_max is the supremum of the minimal range of the sum S_k.
# For N steps, the worst-case minimal range L_max = 2 * (1 - 1/N).
# So, we need 2 * (1 - 1/N) <= 2 * (1 - 10**(-M))
# This simplifies to 1/N >= 10**(-M)
# which is equivalent to 10**M >= N
# M >= log10(N)

# Calculate the minimum required value for M
M_float = math.log10(N)

# Since M must be a positive integer, we take the ceiling of the result.
M = math.ceil(M_float)

print(f"Let N be the number of real numbers a_i, which is {N}.")
print(f"The problem requires finding the smallest integer M such that for any sequence a_i, a valid sequence x_i can be found.")
print(f"This reduces to the condition max_k(S_k) - min_k(S_k) <= 2 * S, where S = 1 - 10**(-M).")
print(f"Let L_max be the supremum over all sequences {a_i} of the minimum possible range of partial sums.")
print(f"The condition becomes L_max <= 2 * S.")
print(f"It's a known result from discrepancy theory that L_max = 2 * (1 - 1/N).")
print(f"So, we must have 2 * (1 - 1/{N}) <= 2 * (1 - 10**(-M)).")
print(f"This simplifies to 1/{N} >= 10**(-M), or 10**M >= {N}.")
print(f"Taking the log base 10 of both sides: M >= log10({N}).")
print(f"log10({N}) = {M_float}")
print(f"Since M must be the smallest positive integer satisfying this, we must have M = ceil({M_float}).")
print(f"The smallest integer M is {M}.")

<<<5>>>