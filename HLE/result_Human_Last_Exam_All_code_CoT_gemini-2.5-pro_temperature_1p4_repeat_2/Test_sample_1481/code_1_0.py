import sympy

# Define the symbols used in the equation
s, s_prime = sympy.symbols("s s'")
kappa_s = sympy.Function('κ')(s)
H_s_s_prime = sympy.Function('H')(s, s_prime)

# The leading order term
term1 = sympy.Abs(s_prime - s)

# The first correction term
term2 = (sympy.Integer(1)/sympy.Integer(24)) * kappa_s**2 * sympy.Abs(s_prime - s)**3

# Higher order terms
big_O = sympy.O((s_prime - s)**4)

# Construct the full asymptotic expression
asymptotic_expansion = sympy.Eq(H_s_s_prime, term1 - term2 + big_O)

# --- Output the result ---
# The task requires printing the equation with each number explicitly.
# We will format the output as a human-readable string.
print("The asymptotic expansion of the generating function H(s, s') is:")
print(f"{H_s_s_prime.name}(s, s') ≈ |s' - s| - (1/24) * κ(s)² * |s' - s|³ + O(|s' - s|⁴)")
