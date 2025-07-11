import math

N = 100000

# The worst-case span W_N is 2 * (1 - 1/N) for even N
W_N = 2 * (1 - 1/N)

# The condition is 2 * (1 - 10**(-M)) >= W_N
# 1 - 10**(-M) >= 1 - 1/N
# 10**(-M) <= 1/N
# -M <= log10(1/N)
# -M <= -log10(N)
# M >= log10(N)

# We are looking for the smallest integer M satisfying this.
M = math.ceil(math.log10(N))

print(f"The number of real numbers is N = {N}")
print(f"The maximum minimal span is W_N = 2 * (1 - 1/N) = {W_N}")
print(f"The condition on M is 2 * (1 - 10**(-M)) >= W_N")
print(f"This simplifies to 10**(-M) <= 1/N, which means M >= log10(N)")
print(f"log10(N) = log10({N}) = {math.log10(N)}")
print(f"The smallest integer M must be {M}.")
