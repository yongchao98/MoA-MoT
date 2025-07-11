import math

# --- FMO Theory Parameters ---
# Energy penalty for an electron-donating methyl group rotating inward instead of outward.
dE_penalty_per_group = 1.6  # kcal/mol

# Number of methyl groups contributing to the energy difference.
num_groups = 2

# Total activation energy difference between the disfavored path (to A) and the favored path (to B).
# Ea(A) - Ea(B)
delta_delta_E = num_groups * dE_penalty_per_group  # kcal/mol

# --- Physical Constants and Conditions ---
R = 1.9872 / 1000  # Gas constant in kcal/mol·K
T = 298.15         # Absolute temperature in Kelvin (25 °C)

# --- Calculation ---
# Calculate the product ratio [A]/[B] using the Boltzmann distribution.
# ratio = exp( - (Ea_A - Ea_B) / RT )
ratio_A_to_B = math.exp(-delta_delta_E / (R * T))
ratio_B_to_A = 1 / ratio_A_to_B

# --- Output the results ---
print("--- FMO Theory Prediction for Product Ratio ---")
print(f"The reaction is a thermal 8-pi electron electrocyclization, which proceeds via a conrotatory mechanism.")
print(f"Torquoselectivity predicts the trans-isomer (B) is favored over the cis-isomer (A).")
print("\n--- Calculation Details ---")
print(f"Energy preference for outward rotation of one methyl group: {dE_penalty_per_group} kcal/mol")
print(f"Total activation energy difference Ea(A) - Ea(B): {num_groups} * {dE_penalty_per_group} = {delta_delta_E:.1f} kcal/mol")
print(f"Temperature (T): {T} K")
print(f"Gas Constant (R): {R:.6f} kcal/mol·K")

# Print the final equation with all numbers plugged in
print("\n--- Final Equation ---")
print(f"Ratio [A]/[B] = exp(-({delta_delta_E:.1f}) / ({R:.6f} * {T}))")
print(f"Ratio [A]/[B] = {ratio_A_to_B:.4f}")
print(f"This means the ratio of trans-isomer B to cis-isomer A is approximately {ratio_B_to_A:.1f} : 1.")

# The problem asks for the ratio of A and B. We will provide the decimal ratio A/B.
final_answer = ratio_A_to_B
print(f"\nThe predicted ratio of A to B is {final_answer:.4f}.")
<<<0.0045>>>