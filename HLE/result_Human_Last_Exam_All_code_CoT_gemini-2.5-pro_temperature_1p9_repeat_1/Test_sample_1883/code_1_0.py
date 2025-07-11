import math

# --- Calculation of the Product Ratio A/B ---

# This script calculates the ratio of cis-isomer A to trans-isomer B
# based on the principles of FMO theory and the Curtin-Hammett principle.

# Define the constants for the calculation.
# ΔΔG‡: The difference in activation free energy between the pathways leading to
# products A and B. This is assumed to be the steric strain of one "inward"
# pointing methyl group in the transition state.
delta_G_diff_kJ_per_mol = 7.0

# R: The ideal gas constant in kJ/(mol·K).
R_kJ_per_mol_K = 0.00831446

# T: The absolute temperature in Kelvin. "Thermal condition" is assumed
# to be 100 degrees Celsius.
T_celsius = 100
T_kelvin = T_celsius + 273.15

# Calculate the term RT, which is the thermal energy.
RT_kJ = R_kJ_per_mol_K * T_kelvin

# Calculate the ratio using the Boltzmann distribution.
# Product A is formed via the lower-energy transition state, so it is the major product.
# Ratio A/B = exp(ΔΔG‡ / RT)
ratio_A_to_B = math.exp(delta_G_diff_kJ_per_mol / RT_kJ)

# --- Output Results ---

print("Prediction of the product ratio for the thermal electrocyclization.")
print("The ratio of A (cis-isomer) to B (trans-isomer) is calculated using the formula: Ratio = exp(ΔΔG‡ / (R * T))\n")

print("Here are the values used in the calculation:")
print(f"Activation Energy Difference (ΔΔG‡): {delta_G_diff_kJ_per_mol} kJ/mol")
print(f"Ideal Gas Constant (R): {R_kJ_per_mol_K:.6f} kJ/(mol*K)")
print(f"Temperature (T): {T_kelvin} K ({T_celsius}°C)")
print("-" * 40)
print("Intermediate Calculation:")
print(f"The thermal energy term (R * T) is {RT_kJ:.4f} kJ/mol.")
print(f"The exponent (ΔΔG‡ / RT) is {delta_G_diff_kJ_per_mol / RT_kJ:.4f}.")
print("-" * 40)
print(f"Final equation with numbers: Ratio = exp({delta_G_diff_kJ_per_mol} / {RT_kJ:.4f})")
print(f"\nThe calculated ratio of product A to product B is {ratio_A_to_B:.2f}.")
print(f"This means for every 1 mole of trans-isomer B, approximately {ratio_A_to_B:.2f} moles of cis-isomer A are produced.")
