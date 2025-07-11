import sympy

# Define the symbolic variable 's' and the function K(s)
s = sympy.Symbol('s')
K = sympy.Function('K')(s)

# Define the polynomials from the controller formula
# Numerator: (s**2 - 1) * K(s) + 4*s**2 + 8*s + 4
# Denominator: (s**2 - 1) - s * K(s)

# We use sympy to pretty-print the expressions.
num_poly_coeff = sympy.poly(4*s**2 + 8*s + 4, s)
den_poly = sympy.poly(s**2 - 1, s)
num_k_coeff = den_poly
den_k_coeff = sympy.poly(-s, s)

print("The set of all proper stabilizing controllers H_2(s) is given by:")
print("H_2(s) = Num(s) / Den(s)")
print("\nwhere K(s) is any stable and proper rational function.")
print("-" * 50)

print("\nNumerator of the controller, Num(s):")
# We explicitly print each term and coefficient as requested.
print(f"Num(s) = ({num_k_coeff.as_expr()}) * K(s) + ({num_poly_coeff.as_expr()})")
print(f"Num(s) = (1*s**2 - 1) * K(s) + (4*s**2 + 8*s + 4)")


print("\nDenominator of the controller, Den(s):")
print(f"Den(s) = ({den_poly.as_expr()}) + ({den_k_coeff.as_expr()}) * K(s)")
print(f"Den(s) = (1*s**2 - 1) - 1*s * K(s)")
print("-" * 50)

# Final Answer for the bot
final_answer_num = "({s_sq_minus_1}) * K(s) + ({four_s_sq})".format(
    s_sq_minus_1="s**2 - 1",
    four_s_sq="4*s**2 + 8*s + 4"
)
final_answer_den = "({s_sq_minus_1}) - s * K(s)".format(
    s_sq_minus_1="s**2 - 1"
)
final_answer = "H_2(s) = ({num}) / ({den})".format(
    num=final_answer_num,
    den=final_answer_den
)
