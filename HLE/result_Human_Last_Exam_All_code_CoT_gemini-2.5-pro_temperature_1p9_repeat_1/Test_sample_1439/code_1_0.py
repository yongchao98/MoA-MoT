import sympy

# Define the coupling constant u and a generic coefficient c1 as symbolic variables
u = sympy.Symbol('u')
c1 = sympy.Symbol('c1')

# The one-loop expression for the anomalous dimension γ_t(u) is linear in u.
# This is the first perturbative correction.
# γ_t(u) ≈ c1*u
gamma_t_one_loop = c1 * u

# The critical exponent ν is related to γ_t(u) by the formula:
# ν(u) = 1 / (2 - γ_t(u))
nu_expression = 1 / (2 - gamma_t_one_loop)

# Perform a Taylor series expansion of ν(u) around u=0 to find the corrections
# to the mean-field value.
nu_series = nu_expression.series(u, 0, 3)

print("The expansion for the critical exponent ν as a function of the coupling u is calculated.")
print(f"The symbolic series is: ν(u) ≈ {nu_series}")
print("\nThis can be broken down into its terms:")

# Zeroth order term (Mean-Field value)
term_0_val = 1/2
print(f"Term 0 (Mean-Field value): {term_0_val}")

# First order term (Initial correction)
# Extract the coefficient for printing
coeff_1 = nu_series.coeff(u, 1)
print(f"Term 1 (First Correction): {coeff_1} * u")

print("\nTo satisfy the output format requirement, the final equation for ν(u) is:")
# The resulting expansion is 1/2 + c1/4 * u + O(u**2). The numbers in the final equation are 1, 2, and 4.
print(f"ν(u) ≈ {1}/{2} + (c1/{4}) * u + ...")

print("\nThe initial non-vanishing contribution to ν beyond its mean-field value is the first-order term.")
print("This contribution is proportional to the coupling constant u raised to the power of 1.")