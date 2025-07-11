import sympy

# Define the symbols for the formula
mu = sympy.Symbol('μ')
e = sympy.Symbol('e')

# Define the indices for the stages involved in the second plateau
# The second plateau is the transition between Stage 2 and Stage 3.
k_final = 2
k_initial = 3

# Construct the chemical potential terms
mu_2 = sympy.Indexed(mu, k_final)
mu_3 = sympy.Indexed(mu, k_initial)

# The formula for the voltage plateau is the difference in chemical potentials
# divided by the elementary charge 'e'.
# V_plateau_2 = (μ_2 - μ_3) / e
# This represents the voltage for the transition from stage 3 to stage 2.
formula = (mu_2 - mu_3) / e

# The question states this second plateau is at approximately 0.13V.
# We will print the full equation as requested.
plateau_voltage = 0.13
print(f"The formula for the second plateau (at {plateau_voltage}V) is:")
print(f"{plateau_voltage} V = ({mu_2} - {mu_3}) / {e}")