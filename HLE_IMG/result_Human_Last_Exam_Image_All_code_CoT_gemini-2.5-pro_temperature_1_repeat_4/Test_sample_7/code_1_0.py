import sympy

# Define the symbols for the equation
V, mu_1, mu_2, mu_Li, e = sympy.symbols('V μ_1 μ_2 μ_Li e')

# Construct the equation based on the derivation
# The second plateau voltage 'V' is approximated by a formula involving
# the chemical potentials of stage 1 (mu_1), stage 2 (mu_2),
# the Li reference (mu_Li), and the elementary charge (e).
# The derivation assumes an averaging model for the plateau potential
# and linear extrapolation for the chemical potentials of the stages.
equation = sympy.Eq(V, (mu_1 - 3*mu_2 + 2*mu_Li) / (2*e))

# Print the formula in a human-readable format
# The 'pretty' print function from sympy can format it nicely,
# but for a standard python print, we'll format it manually.
# We explicitly state each number in the formula as requested.
# V = (1*μ_1 - 3*μ_2 + 2*μ_Li) / (2*e)
print("The formula for the second plateau is:")
print("V = (μ_1 - 3*μ_2 + 2*μ_Li) / (2*e)")