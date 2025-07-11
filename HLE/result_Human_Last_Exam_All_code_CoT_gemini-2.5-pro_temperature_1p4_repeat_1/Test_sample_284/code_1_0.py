# The dimension of the space is n.
n = 3

# The problem is about the Fourier extension operator for the moment curve
# gamma(t) = (t, t^2, ..., t^n) in R^n.
# The question asks for the critical exponent p_c such that for p < p_c,
# no non-zero L^p function has Fourier support on the curve, while for p >= p_c,
# such functions can exist.
# This critical exponent is given by the formula p_c = n(n+1)/(n-1).

# Calculate the numerator and the denominator of the formula.
numerator = n * (n + 1)
denominator = n - 1

# Calculate the critical exponent p_c.
p_c = numerator / denominator

# The question asks for the largest possible value of p for which the statement holds.
# This value is the supremum of the set {p | statement holds}, which is p_c.

# Print the calculation step-by-step as requested.
print(f"The dimension is n = {n}")
print(f"The formula for the critical exponent p_c is n*(n+1)/(n-1)")
print(f"Numerator: n * (n + 1) = {n} * ({n} + 1) = {numerator}")
print(f"Denominator: n - 1 = {n} - 1 = {denominator}")
print(f"The largest possible value of p is p_c = {numerator} / {denominator} = {p_c}")