import math

# Number of steps
N = 100000

# The problem is to find the smallest integer M such that for any sequence a_i,
# we can find x_i satisfying the conditions.
# This boils down to the inequality W_max <= 2 - 2 * 10**(-M),
# where W_max is the supremum of the minimal range of partial sums.
# For large N, W_max can be approximated by 2 - 1/N.

# So we have the inequality: 2 - 1/N <= 2 - 2 * 10**(-M)
# Which simplifies to: 2 * 10**(-M) <= 1/N
# 10**(-M) <= 1 / (2 * N)
# Taking log10 on both sides: -M <= log10(1 / (2 * N))
# -M <= -log10(2 * N)
# M >= log10(2 * N)
# M >= log10(2) + log10(N)

log10_N = math.log10(N)
log10_2 = math.log10(2)

M_min = log10_2 + log10_N
M = math.ceil(M_min)

print(f"Let N = {N}.")
print("The maximum minimal range, W_max, is approximated by 2 - 1/N for large N.")
print(f"W_max ≈ 2 - 1/{N} = {2 - 1/N}")
print("We must satisfy the condition: W_max <= 2 - 2 * 10**(-M)")
print(f"So, 2 - 1/{N} <= 2 - 2 * 10**(-M)")
print(f"2 * 10**(-M) <= 1/{N}")
print(f"10**(-M) <= 1 / (2 * {N})")
print(f"-M <= log10(1 / (2 * {N}))")
print(f"M >= log10(2 * {N})")
print(f"M >= log10(2) + log10({N})")
print(f"M >= {log10_2} + {log10_N}")
print(f"M >= {M_min}")
print(f"Since M must be the smallest positive integer, M = {M}")