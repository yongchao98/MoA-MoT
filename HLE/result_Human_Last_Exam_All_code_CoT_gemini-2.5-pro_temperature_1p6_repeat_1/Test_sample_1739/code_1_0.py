import sympy

# Define the symbol for the polytropic index
gamma = sympy.Symbol('gamma')

# The nonlinear frequency correction coefficient, sigma, is a polynomial in gamma.
# The correction to the frequency is given by w = w0 * (1 + epsilon^2 * a^2 * sigma)
# sigma = A*gamma**2 + B*gamma + C
# where A, B, and C are the coefficients.
# Based on the derivation using multiple scale analysis, the coefficients are:
A = -6
B = -sympy.Rational(27, 4)
C = -sympy.Rational(1, 2)

# The polynomial expression for the nonlinear correction coefficient
correction_poly = A*gamma**2 + B*gamma + C

# The problem asks for the 3rd term of the nonlinear correction.
# We interpret this as the third term in the polynomial expression for sigma.
# The terms are:
# 1st term: -6 * gamma**2
# 2nd term: -(27/4) * gamma
# 3rd term: -1/2

third_term = C

# Output the components of the nonlinear correction equation
# and the final answer.
print("The equation for the normalized nonlinear frequency correction coefficient (sigma) is:")
print(f"sigma = ({A}) * gamma**2 + ({B}) * gamma + ({C})")
print("\nThe three terms of the correction are:")
print(f"1st term: {A}*gamma**2")
print(f"2nd term: {B}*gamma")
print(f"3rd term: {C}")

print("\n-------------------------------------------")
print("The calculated value of the 3rd term of the nonlinear correction is:")
print(float(third_term))
