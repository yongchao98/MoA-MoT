import sympy

# This script calculates the second-order energy shift due to the relativistic
# kinetic energy correction for a hydrogen atom in the state n=3, l=2.

# --- Symbolic setup ---
# Define fundamental constants as symbols for algebraic manipulation.
m_e, c, alpha = sympy.symbols('m_e c alpha')

# Define the quantum numbers for the state in question.
n_val = 3
l_val = 2

# --- Derivations and Calculations ---

# The direct calculation of the second-order energy shift involves an infinite sum over
# all other states and is extremely complex. A standard method for this problem is
# to use the Unsold approximation.
#
# Plan:
# 1. Calculate the first-order energy shift, ΔE^(1).
# 2. Use the approximation ΔE^(2) ≈ -(ΔE^(1))^2 / (2|E_n|) to find the result.

# Step 1: Calculate the first-order energy shift ΔE^(1).
# The formula for the first-order shift for this perturbation is:
# ΔE^(1)_nl = - (E_n^2 / (2*m_e*c^2)) * (4n/(l + 1/2) - 3)
# This is derived from the expectation value <H'>.
n, l = sympy.symbols('n l')
E_n_formula = -(m_e * c**2 * alpha**2) / (2 * n**2)
delta_E_1_formula = - (E_n_formula**2 / (2*m_e*c**2)) * (4*n/(l + sympy.Rational(1, 2)) - 3)

# Substitute n=3, l=2 into the formula.
delta_E_1 = sympy.simplify(delta_E_1_formula.subs([(n, n_val), (l, l_val)]))
# The symbolic result is ΔE^(1) = (-1/360) * m_e * c^2 * alpha^4

# Step 2: Calculate the second-order energy shift ΔE^(2).
# The energy of the unperturbed state n=3 is E_3.
E_n_val = E_n_formula.subs(n, n_val)
abs_E_n_val = sympy.Abs(E_n_val) # |E_3| = m_e * c^2 * alpha^2 / 18

# Apply the Unsold approximation formula.
delta_E_2 = sympy.simplify(- (delta_E_1**2) / (2 * abs_E_n_val))

# --- Output the Final Result ---
# The result of the symbolic calculation is ΔE^(2) = (-1/14400) * m_e * c^2 * alpha^6.
# We extract the coefficient for clear printing.
coeff = delta_E_2 / (m_e * c**2 * alpha**6)
num, den = coeff.as_numer_denom()

print("The second-order relativistic energy shift ΔE^(2) for the state n=3, l=2 is calculated as follows:")
print("The final equation for the energy shift is:")
print(f"ΔE^(2) = ({num} / {den}) * m_e * c^2 * α^6")