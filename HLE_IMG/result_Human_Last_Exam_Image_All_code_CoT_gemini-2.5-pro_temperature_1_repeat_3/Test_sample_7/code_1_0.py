# The goal is to derive and present the formula for the second voltage plateau
# of a graphite anode based on the provided information.

# The analysis leads to the conclusion that the second plateau, which corresponds
# to the formation of the Stage 2 Li-Graphite intercalated compound, is
# governed by the chemical potential of that stage, denoted as μ_2.

# The fundamental relationship between the half-cell voltage (V),
# the chemical potential (μ) of the intercalating species, and its charge (e)
# is V = -μ / e.

# Therefore, the formula for the second plateau is V = -μ_2 / e.

# The code below prints this final formula and its constituent parts.

# Define the symbols for the components of the equation
sign = "-"
numerator = "μ_2"
denominator = "e"
variable = "V"

# Print the final equation
print(f"The simple formula that best approximates the second plateau is:")
print(f"{variable} = {sign}{numerator} / {denominator}")

# As requested, output each component of the final equation
print("\nComponents of the final equation:")
print(f"'{variable}': Represents the voltage of the second plateau (approx. 0.13V).")
print(f"'{sign}': The negative sign, which comes from the thermodynamic relationship between Gibbs free energy and electrode potential.")
print(f"'{numerator}': Represents the chemical potential of intercalating a single Li+ ion at (almost full) stage 2.")
print(f"'{denominator}': Represents the elementary charge of the Li+ ion.")