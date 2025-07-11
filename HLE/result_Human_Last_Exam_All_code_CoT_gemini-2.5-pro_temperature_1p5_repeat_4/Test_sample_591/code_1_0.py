import sys

# Define the components of the formula as a dictionary for clarity.
# Using unicode characters for Greek letters to match the prompt.
variables = {
    "κ (kappa)": "A dimensionless parameter representing a ratio of thresholds",
    "ρ (rho)": "An offset constant in the Hebbian equation",
    "τ_u (tau_u)": "Time constant for the postsynaptic accumulator",
    "τ_v (tau_v)": "Time constant for the presynaptic accumulator",
    "φ (phi)": "A scaling constant for the presynaptic accumulator"
}

# The formula for kappa as a string.
# Note: In Python code, you would use variable names like 'rho', 'tau_u', etc.
# Here, we print the symbolic representation.
# To show the "final equation" and its "numbers" (in this case, variables),
# we construct the string from its components.
rho = "ρ"
tau_u = "τ_u"
tau_v = "τ_v"
phi = "φ"

numerator = f"-{rho} * {tau_u} * ({tau_u} + {tau_v})"
denominator = f"{phi} * {tau_v}"
kappa_formula = f"κ = ({numerator}) / ({denominator})"

# Print the final result
print("The definition of the parameter κ in the given model is:")
print(kappa_formula)

print("\nWhere the variables in the equation are:")
for var, description in variables.items():
    print(f"- {var}: {description}")

# Provide the qualitative interpretation from the source paper
explanation = (
    "In the context of the model, κ represents the ratio of the depression "
    "threshold for postsynaptic activity (set by the ρ parameter) to the mean "
    "potential induced by presynaptic activity."
)
print(f"\nInterpretation: {explanation}")