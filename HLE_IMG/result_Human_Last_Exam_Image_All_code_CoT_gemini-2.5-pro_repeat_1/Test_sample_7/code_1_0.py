import sympy

# Define the symbols for the chemical potentials and the elementary charge
mu_1 = sympy.Symbol('μ_1')
mu_2 = sympy.Symbol('μ_2')
e = sympy.Symbol('e')
V = sympy.Symbol('V')

# The problem describes two main voltage plateaus for Li intercalation in graphite.
# The first plateau (approx. 0.09V) corresponds to the Stage 2 -> Stage 1 transition.
# The second plateau (approx. 0.13V) corresponds to the transition from a dilute phase to Stage 2.

# A simple approximation relates the plateau voltage to the chemical potential of the final stage,
# i.e., V_plateau2 ≈ -μ_2 / e. However, the prompt asks for a formula using different μ_k.

# A more advanced model from the physics of staging compounds approximates the
# chemical potential change for the second plateau (dilute -> stage 2) as Δμ ≈ 2*μ_2 - μ_1.
# The voltage is given by V = -Δμ / e.

# This results in the formula V ≈ -(2*μ_2 - μ_1) / e.
# We will print this formula. The numbers '1' and '2' are part of the symbols μ_1 and μ_2.

# Construct the equation
equation = sympy.Eq(V, (mu_1 - 2 * mu_2) / e)

# Print the formula in a readable format
# The instruction "output each number in the final equation" refers to the coefficients and indices.
print("The formula for the second plateau is:")
print(f"V = (μ_1 - 2 * μ_2) / e")

# This can also be written as:
print("\nOr, alternatively:")
print(f"V = -(2 * μ_2 - μ_1) / e")
