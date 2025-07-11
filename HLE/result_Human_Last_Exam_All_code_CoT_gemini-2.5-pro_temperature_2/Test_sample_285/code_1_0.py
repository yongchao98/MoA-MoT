# The dimension of the parameter space for the integral I(a_1, ..., a_9)
N = 9

# The slowest decay exponent sigma_min for the integral I(a).
# This is determined by the highest order terms in the phase polynomial,
# which are cubic (degree 3).
# For a polynomial of degree d in n variables, the decay exponent is typically n/d.
# However, for an integral on a domain with boundaries, we must consider
# contributions from boundary degeneracies.
# The analysis shows the slowest decay comes from terms like a_6*x^3, where the
# phase has a highly degenerate critical point at the boundary, giving a
# decay exponent of 1/3.
sigma_min = 1/3

# The Lp norm of I(a) is finite if the integral of |I(a)|^p converges.
# In hyperspherical coordinates, the integral is proportional to:
# integral from 0 to inf of r^(N-1) * |I(r_vec)|^p dr
# The integrand behaves as r^(N-1) * (r^(-sigma_min))^p = r^(N-1 - p*sigma_min).
# This integral diverges if the exponent is >= -1.
# N - 1 - p * sigma_min >= -1
# which simplifies to p <= N / sigma_min.
p_max = N / sigma_min

print(f"The dimension of the parameter space is N = {N}.")
print(f"The slowest decay exponent of the integral is sigma_min = 1/{int(1/sigma_min)}.")
print("The condition for the integral not to be in L^p is given by the inequality:")
print(f"p * sigma_min <= N")
print(f"p * (1/{int(1/sigma_min)}) <= {N}")
print(f"p <= {N} * {int(1/sigma_min)}")
print(f"p <= {p_max}")
print("\nThe largest value of p for which the function I is not in L^p(R^9) is the maximum value satisfying this inequality.")
print(f"The largest p is {p_max}.")