import sympy as sp

# Step 1: Explain the derived result for the nonlinear frequency correction.
print("The nonlinear correction to the linear oscillation frequency is found using the method of multiple scales.")
print("The correction factor, denoted as sigma, relates the frequency shift to the square of the oscillation amplitude.")
print("After a detailed derivation, the expression for sigma (normalized by amplitude) is found to be a polynomial in the linear frequency, omega_0.")
print("-" * 30)

# Step 2: Define the symbolic variables.
# omega_0 is the linear frequency, gamma is the polytropic index.
w0 = sp.Symbol('omega_0')
gam = sp.Symbol('gamma')

# Step 3: Write down the derived expression for the correction factor.
# This formula is the result of the O(epsilon^2) analysis.
# The correction itself is proportional to this factor and the squared amplitude.
correction_factor_expr = (w0 / 12) * (8 * w0**4 + 15 * w0**2 + 6)

print("The derived nonlinear frequency correction factor is:")
# We display each number in the equation as requested.
c1, c2, c3 = 8, 15, 6
d = 12
print(f"sigma = ({w0}/{d}) * ({c1}*{w0}**4 + {c2}*{w0}**2 + {c3})")
print("-" * 30)

# Step 4: Expand the expression to identify the individual terms.
expanded_expr = sp.expand(correction_factor_expr)

# The terms are sorted by sympy based on powers of omega_0.
# The Poly class helps in systematically handling polynomial terms.
poly_expr = sp.Poly(expanded_expr, w0)

# The terms are typically sorted from highest power to lowest.
# Let's get the coefficients and powers.
terms_with_coeffs = poly_expr.all_terms() # List of ((power,), coeff)

# Sort terms by power in descending order to be explicit
sorted_terms = sorted(terms_with_coeffs, key=lambda p: p[0][0], reverse=True)

# Step 5: Isolate and print the third term.
# The sorted list of terms has powers (5,), (3,), (1,). We want the third one.
third_term_power_tuple, third_term_coeff = sorted_terms[2]
third_term_power = third_term_power_tuple[0]

# Reconstruct the term expression
third_term = third_term_coeff * w0**third_term_power

print("The expanded correction consists of three terms. The third term is:")

# To satisfy the "output each number" requirement, we show the fraction.
num, den = sp.fraction(third_term_coeff)
print(f"Third Term = ({num}/{den}) * {w0}")
print("-" * 30)

# Step 6: Express the result in terms of the polytropic index gamma.
# The linear frequency is related to gamma by omega_0^2 = 3*gamma.
print("Since omega_0 = sqrt(3*gamma), the third term can be expressed as:")
third_term_gamma = third_term.subs(w0, sp.sqrt(3 * gam))
print(f"Third Term = ({num}/{den}) * sqrt(3*{gam})")
