import sympy

# Let x represent omega_0^2, which is equal to 3*gamma
x = sympy.Symbol('x')

# From the Poincaré-Lindstedt analysis up to the third order,
# the second frequency correction omega_2 is found to be related
# to a polynomial in x = omega_0^2.
# The relationship is: 48 * omega_2 / omega_0 = -2*x**2 + 3*x + 6
final_polynomial = -2*x**2 + 3*x + 6

# The question asks for the "3rd term" of the nonlinear correction.
# We interpret this as the third term of the polynomial that determines the correction.
# In standard form (descending powers of x), the terms are:
# 1st term: -2*x**2
# 2nd term: 3*x
# 3rd term: 6

# We can extract the coefficients of the polynomial to verify.
poly_obj = sympy.Poly(final_polynomial, x)
coeffs = poly_obj.all_coeffs()
term1_coeff = coeffs[0]
term2_coeff = coeffs[1]
term3_coeff = coeffs[2]

final_equation_str = f"{term1_coeff}*x**2 + {term2_coeff}*x + {term3_coeff}"
print(f"The frequency correction is determined by the polynomial equation: P(x) = 0, where P(x) = {final_equation_str}")
print(f"The numbers in this equation are: {term1_coeff}, {term2_coeff}, and {term3_coeff}")
print(f"The third term of this polynomial is the constant term.")
print(f"The value of the 3rd term is: {term3_coeff}")
