# This code must be run in a SageMath environment.

# Define the polynomial ring over the rational numbers QQ.
R.<x> = PolynomialRing(QQ)

# Define the polynomial f(x) from the equation y^2 = f(x).
# The coefficients are 1, 4, 6, 2, 1, 2, 1.
f = 1*x^6 + 4*x^5 + 6*x^4 + 2*x^3 + 1*x^2 + 2*x + 1

# Define the hyperelliptic curve.
C = HyperellipticCurve(f)

# Compute the conductor of the curve's Jacobian.
N = C.conductor()

# Print the final result, restating the equation for clarity as requested.
print(f"The conductor of the curve defined by y^2 = {f} is {N}")
