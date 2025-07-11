# This script prints the definition of the parameter kappa (κ).
# The equation for κ is derived from the steady-state analysis of the provided plasticity model.

print("The definition of κ is given by the following equation:")

# The instruction asks to output each "number" (symbol) in the final equation.
# The following lines print the formula piece by piece to meet this requirement.
print("κ = -(", end="")
print("μ", end="")
print(" + ", end="")
print("ρ", end="")
print("/", end="")
print("φ", end="")
print(") * (", end="")
print("τ_u", end="")
print(" + ", end="")
print("τ_v", end="")
print(")")