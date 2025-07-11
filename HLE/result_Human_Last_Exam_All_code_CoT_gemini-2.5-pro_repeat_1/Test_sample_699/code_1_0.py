import math

# --- Given Parameters and Assumptions ---
b = -0.09  # Basquin exponent

# Cycle fractions for each stress level
f1, f2, f3 = 0.70, 0.20, 0.10

# Stress levels are sigma_1 = sigma_e, sigma_2 = 1.1*sigma_e, sigma_3 = 1.2*sigma_e
stress_ratio_1 = 1.0
stress_ratio_2 = 1.1
stress_ratio_3 = 1.2

# Assumption: The problem is unsolvable without knowing the life at the endurance limit (N1).
# We assume a standard engineering value for the life at the endurance limit.
N1 = 1_000_000  # Assumed life in cycles at sigma_e

# --- Calculations ---

# Step 1: Calculate lives N2 and N3 based on N1 using the Basquin relationship.
# N_i = N_1 * (sigma_i / sigma_1)^(1/b)
# Note: Since b is negative, a higher stress results in a lower life.
N2 = N1 * (stress_ratio_2)**(1/b)
N3 = N1 * (stress_ratio_3)**(1/b)

# Step 2: Apply the Palmgren-Miner rule to find the total fatigue life, Nf.
# 1/Nf = f1/N1 + f2/N2 + f3/N3
total_damage_per_cycle = (f1 / N1) + (f2 / N2) + (f3 / N3)
Nf = 1 / total_damage_per_cycle

# --- Output Results ---
print("Fatigue Life Estimation based on Palmgren-Miner Rule")
print("-" * 55)
print(f"Assumption: Life at the endurance limit (N1) = {N1:,.0f} cycles.")
print(f"Basquin exponent (b) = {b}")
print("-" * 55)

# Per the user request, output each number in the final equation.
# The final equation is: 1/Nf = f1/N1 + f2/N2 + f3/N3
print("The fatigue life equation is: 1/Nf = f1/N1 + f2/N2 + f3/N3")
print("\nWhere the values are:")
print(f"f1 = {f1:.2f}, N1 = {N1:15,.3f} cycles")
print(f"f2 = {f2:.2f}, N2 = {N2:15,.3f} cycles")
print(f"f3 = {f3:.2f}, N3 = {N3:15,.3f} cycles")
print("-" * 55)

# Final calculated life
print(f"Calculated Total Fatigue Life (Nf) = {Nf:,.3f} cycles")

# Format the final answer as Aeb with 3 decimal places
formatted_Nf = f"{Nf:.3e}".replace('e+0', 'e')
print(f"\nFinal Answer (in Aeb format): {formatted_Nf}")

# The final answer is enclosed below as requested
# <<<4.966e5>>>