import math

# The problem is to find the smallest positive integer M such that for any
# a_1, ..., a_N in [0,1], there exist x_0, ..., x_N in [-S, S]
# with S = 1 - 10^(-M) and |x_{i-1} - x_i| = a_i.
# This can be shown to be equivalent to the condition that the maximum possible value
# of the minimal range of partial sums is less than or equal to 2S.
# Let C(N) be this maximum minimal range. The condition is C(N) <= 2*S.
# For large N, a very good approximation for C(N) is 2 - 2 / sqrt(pi * N / 2).
# So, we have the inequality: 2 - 2/sqrt(pi*N/2) <= 2 * (1 - 10^(-M))

# Let's solve for M.
# 1 - 1/sqrt(pi*N/2) <= 1 - 10^(-M)
# 10^(-M) <= 1/sqrt(pi*N/2)
# -M <= log10(1/sqrt(pi*N/2))
# -M <= -log10(sqrt(pi*N/2))
# M >= log10(sqrt(pi*N/2))

N = 100000

# Calculate the value of log10(sqrt(pi*N/2))
val = math.log10(math.sqrt(math.pi * N / 2.0))

# The smallest integer M must be ceiling(val)
M = math.ceil(val)

print(f"The number of real numbers is N = {N}")
print(f"The condition is R_max <= 2 * (1 - 10^(-M))")
print(f"The maximum minimal range R_max for N={N} is approximately C({N})")
print(f"C({N}) approx 2 - 2/sqrt(pi*N/2)")
print(f"This leads to the inequality M >= log10(sqrt(pi*N/2))")
print(f"Calculating the right hand side:")
print(f"pi * N / 2 = {math.pi * N / 2.0}")
print(f"sqrt(pi * N / 2) = {math.sqrt(math.pi * N / 2.0)}")
print(f"log10(sqrt(pi * N / 2)) = {val}")
print(f"The smallest integer M satisfying M >= {val} is {M}.")
print(f"The final equation is M >= log10(sqrt(pi*{N}/2))")
print(f"M >= log10(sqrt({math.pi * N / 2.0}))")
print(f"M >= {val}")
print(f"So the smallest integer M is {M}")
