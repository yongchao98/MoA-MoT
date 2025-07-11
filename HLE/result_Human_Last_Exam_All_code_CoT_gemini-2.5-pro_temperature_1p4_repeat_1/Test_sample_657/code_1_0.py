import sympy

# Setup symbolic variables
s = sympy.Symbol('s')
# Define K(s) as a symbolic function
K_func = sympy.Function('K')
K = K_func(s)

# The initial transfer function H_1(s) = s / (s^2 - 1)
H1 = s / (s**2 - 1)

# Step 1: Find a stable coprime factorization H_1(s) = N(s)/D(s).
# The plant H_1(s) has an unstable pole at s=1. We must choose D(s) to have a zero at s=1.
# A standard choice is a Blaschke product (an all-pass filter).
# We choose D(s) = (s-1)/(s+1). It is stable, proper, and has a zero at s=1.
D = (s - 1) / (s + 1)
# Now, we find N(s) from H_1 = N/D, so N = H_1 * D.
N = sympy.cancel(H1 * D)

# Step 2: Solve the Bezout identity X(s)N(s) + Y(s)D(s) = 1.
# From algebraic derivation, a particular stable, proper solution is:
X0 = 4
Y0 = (s - 1) / (s + 1)

# We can verify this solution:
# sympy.simplify(X0 * N + Y0 * D) should yield 1.

# Step 3 & 4: Construct the controller using the Youla-Kucera formula and simplify.
# H_2(s) = (X_0(s) + D(s)K(s)) / (Y_0(s) - N(s)K(s))
H2_numerator_frac = X0 + D * K
H2_denominator_frac = Y0 - N * K

# To simplify the expression, we can multiply the numerator and denominator by (s+1)**2,
# which is the denominator of N(s) and the square of the denominator of D(s).
common_multiplier = (s + 1)**2

# Perform the multiplication and expand the polynomials
numerator_poly = sympy.expand(H2_numerator_frac * common_multiplier)
denominator_poly = sympy.expand(H2_denominator_frac * common_multiplier)

# The final controller is of the form: (A(s)*K(s) + B(s)) / (C(s)*K(s) + D_term(s))
# We extract each polynomial part.
num_k_coeff = sympy.poly(numerator_poly.coeff(K), s)
num_const_term = sympy.poly(numerator_poly.subs(K, 0), s)

den_k_coeff = sympy.poly(denominator_poly.coeff(K), s)
den_const_term = sympy.poly(denominator_poly.subs(K, 0), s)

# Print the final result in a structured way.
print("The set of all proper stabilizing controllers H_2(s) is given by:")

final_eq_str = f"H_2(s) = (({num_k_coeff.as_expr()}) * K(s) + ({num_const_term.as_expr()})) / (({den_k_coeff.as_expr()}) * K(s) + ({den_const_term.as_expr()}))"
print(final_eq_str)
print("\nWhere K(s) is any stable and proper transfer function.")

# As requested, outputting each number in the final equation.
# This is done by listing the coefficients of each polynomial part of the equation.
print("\nCoefficients of the polynomials in the equation (for descending powers of s):")

A_coeffs = num_k_coeff.all_coeffs()
B_coeffs = num_const_term.all_coeffs()
C_coeffs = den_k_coeff.all_coeffs()
D_coeffs = den_const_term.all_coeffs()

print(f"Polynomial multiplying K(s) in Numerator: '{num_k_coeff.as_expr()}' -> coefficients are [{', '.join(map(str, A_coeffs))}]")
print(f"Constant polynomial in Numerator: '{num_const_term.as_expr()}' -> coefficients are [{', '.join(map(str, B_coeffs))}]")
print(f"Polynomial multiplying K(s) in Denominator: '{den_k_coeff.as_expr()}' -> coefficients are [{', '.join(map(str, C_coeffs))}]")
print(f"Constant polynomial in Denominator: '{den_const_term.as_expr()}' -> coefficients are [{', '.join(map(str, D_coeffs))}]")