# The problem asks for the largest p such that the integral of |I|^p diverges.
# Our derivation leads to a condition for divergence based on an exponent.

# d_codim: codimension of the singularity manifold.
# For a locally cubic singularity in 2D, we need to cancel 5 lower-order terms (2 linear, 3 quadratic)
# in the Taylor series of the phase. This imposes 5 constraints on the 9 coefficients.
d_codim = 5

# n: dimension of the ambient space of coefficients.
n = 9

# d_sing: dimension of the singularity manifold.
d_sing = n - d_codim

# decay_exp: the decay exponent of the integral for coefficients on this manifold.
# For a generic cubic phase in 2 variables, the decay is lambda^(-2/3).
decay_exp = 2.0 / 3.0

# The L^p norm integral over the singularity manifold is approximately
# integral(R^(-decay_exp * p) * R^(d_sing - 1) dR)
# where R is the magnitude of the leading coefficients.
# The exponent of R in the integrand is (d_sing - 1 - decay_exp * p).
# The integral diverges if this exponent is >= -1.
# d_sing - 1 - decay_exp * p >= -1
# d_sing >= decay_exp * p
# p <= d_sing / decay_exp

p_critical = d_sing / decay_exp

print(f"The dimension of the ambient coefficient space is n = {n}.")
print(f"The number of constraints to create a generic cubic singularity at a point is 5.")
print(f"The dimension of the singularity manifold is d_sing = n - 5 = {d_sing}.")
print(f"The decay exponent of the integral for this singularity is decay_exp = {decay_exp:.4f}.")
print(f"The critical exponent p is given by the condition: p <= d_sing / decay_exp.")
print(f"So, p <= {d_sing} / ({decay_exp:.4f})")
print(f"p <= {p_critical}")
print("The largest p for which I is not in L^p is 6.")

# Final answer calculation
final_answer = int(p_critical)