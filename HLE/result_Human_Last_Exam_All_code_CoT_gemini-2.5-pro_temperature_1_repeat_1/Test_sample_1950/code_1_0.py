import sympy

# Define symbols for the fidelities
F1, F2 = sympy.symbols('F1 F2')

# Coefficients from the problem description
A_coeff = (8*F1 - 1) * (4*F2 - 1)
B_coeff = (8*F1 - 1) * (1 - F2)
C_coeff = (1 - F1) * (4*F2 - 1)
D_coeff = (1 - F1) * (1 - F2)

# Contribution factors calculated from the protocol analysis
X1_factor = 1
X2_factor = 1
X3_factor = 1
X4_factor = 2

# Total product of fidelity and success probability
total_product = A_coeff * X1_factor + B_coeff * X2_factor + C_coeff * X3_factor + D_coeff * X4_factor

# The expression is divided by (7 * 3) = 21
final_expression = sympy.simplify(total_product / 21)

# Extract coefficients for printing
poly_expr = sympy.Poly(final_expression, F1, F2)
c_f1f2 = poly_expr.coeff_monomial(F1*F2)
c_f1 = poly_expr.coeff_monomial(F1)
c_f2 = poly_expr.coeff_monomial(F2)
c_const = poly_expr.coeff_monomial(1)
denominator = sympy.fraction(final_expression)[1]


print("The product of the successful output fidelity and the success probability is:")
print(f"({c_f1f2}*F1*F2 + {c_f1}*F1 + {c_f2}*F2 + {c_const}) / {denominator}")