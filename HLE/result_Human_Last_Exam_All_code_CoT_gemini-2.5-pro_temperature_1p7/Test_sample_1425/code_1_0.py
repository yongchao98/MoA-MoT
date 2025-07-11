import math

# This script will derive and print the formula for the grand canonical partition function Z
# using the path integral formalism.

# --- Define the symbols used in the formula as string variables for clarity ---

# Z: The Grand Canonical Partition Function
Z = "Z"

# ∫ D[φ]: The path integral over all field configurations φ
path_integral_symbol = "∫ D[φ]"

# exp(...): The exponential function
exponential_function = "exp"

# S_E[φ]: The Euclidean action for the system with Hamiltonian H
euclidean_action = "S_E[φ]"

# μ: The chemical potential
chemical_potential = "μ"

# ∫[0 to β] dτ: Integral over imaginary time τ from 0 to β
imaginary_time_integral = f"∫[0 to {chr(946)}] dτ" # chr(946) is the unicode for beta (β)

# N[φ]: The number operator as a functional of the fields φ
number_operator = "N[φ]"

# --- Construct the terms inside the exponential ---
# The term originating from the Euclidean action
action_term = f"-{euclidean_action}"

# The term originating from the chemical potential coupling to the number operator
# This includes the numbers 0 and β as the integration limits.
chemical_potential_term = f"{chemical_potential} {imaginary_time_integral} {number_operator}"

# Combine the terms for the full expression in the exponent
exponent_expression = f"({action_term} + {chemical_potential_term})"

# --- Assemble the final equation ---
final_equation = f"{Z} = {path_integral_symbol} {exponential_function}{exponent_expression}"

# --- Print the explanation and the final formula ---
print("The grand canonical partition function (Z) is derived from the canonical one by replacing the Hamiltonian H with H' = H - μN.")
print("This modifies the Euclidean action S_E in the path integral exponent.")
print("\nThe original Euclidean action is:")
print(f"  {euclidean_action} = ∫[0 to {chr(946)}] dτ L_E   (where L_E is the Euclidean Lagrangian)")
print("\nThe new effective action becomes:")
print(f"  S'_E = {euclidean_action} - {chemical_potential} {imaginary_time_integral} {number_operator}")
print("\nSubstituting this into the path integral Z = ∫ D[φ] exp(-S'_E) gives the final result:")
print("-" * 60)
print("Final Formula:")
print(final_equation)
print("-" * 60)
