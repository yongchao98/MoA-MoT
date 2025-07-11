import math

# Step 1: Define given parameters and assumptions
b = -0.09  # Basquin exponent
sigma_u_over_sigma_e = 2.0  # Assumed ratio of ultimate strength to endurance limit

# Loading conditions. Note: The load at sigma_e is ignored as it causes no damage.
damaging_loads = [
    {'life_fraction': 0.20, 'stress_factor': 1.1},  # 20% of life at 1.1 * sigma_e
    {'life_fraction': 0.10, 'stress_factor': 1.2}   # 10% of life at 1.2 * sigma_e
]

# Step 2 & 3: Calculate N_e, the theoretical cycles to failure at the endurance limit.
# The S-N curve is defined by two points: (sigma_u, 0.5) and (sigma_e, N_e).
# From sigma_u = C * (0.5)^b and sigma_e = C * N_e^b, we can derive:
# N_e = 0.5 * (sigma_u / sigma_e)^(-1/b)
inv_b = 1.0 / b
N_e = 0.5 * (sigma_u_over_sigma_e)**(-inv_b)

# Step 4, 5, 6: Apply Palmgren-Miner rule and solve for N_total
# The rule is Sum(n_i / N_fi) = 1, which becomes N_total * Sum(p_i / N_fi) = 1.
# We calculate the denominator term: Sum(p_i / N_fi)
denominator = 0.0

# Calculate N_fi for each damaging load and the corresponding damage term
N_f_values = []
damage_terms = []
for load in damaging_loads:
    p_i = load['life_fraction']
    f_i = load['stress_factor']
    
    # Calculate N_fi for this stress level: N_fi = N_e * (stress_factor)^(1/b)
    N_fi = N_e * (f_i)**(inv_b)
    N_f_values.append(N_fi)
    
    # Calculate the damage contribution for the denominator
    damage_contribution = p_i / N_fi
    damage_terms.append(damage_contribution)
    denominator += damage_contribution

# Calculate the total fatigue life
N_total = 1.0 / denominator

# Step 7: Print the results and the final equation
print("--- Calculation Details ---")
print(f"Basquin exponent (b): {b}")
print(f"Assumed sigma_u / sigma_e ratio: {sigma_u_over_sigma_e}")
print(f"Theoretical cycles to failure at endurance limit (N_e): {N_e:.2f}")
print("-" * 30)
print("Palmgren-Miner Damage Calculation (ignoring non-damaging cycles at sigma_e):")
print(f"Cycles to failure at 1.1*sigma_e (N_f1): {N_f_values[0]:.2f}")
print(f"Cycles to failure at 1.2*sigma_e (N_f2): {N_f_values[1]:.2f}")
print("-" * 30)
print("Final Equation for Total Life (N_total):")
print(f"N_total = 1 / ( (Life Fraction 1 / N_f1) + (Life Fraction 2 / N_f2) )")
print(f"N_total = 1 / ( ({damaging_loads[0]['life_fraction']} / {N_f_values[0]:.2f}) + ({damaging_loads[1]['life_fraction']} / {N_f_values[1]:.2f}) )")
print(f"N_total = 1 / ( {damage_terms[0]:.4e} + {damage_terms[1]:.4e} )")
print(f"N_total = 1 / ( {denominator:.4e} )")
print("-" * 30)
print(f"Estimated Fatigue Life (N_total): {N_total:.3f} cycles")

# Format the final answer as requested
formatted_output = f"{N_total:.3e}"
print(f"\nFinal Answer in Aeb format: {formatted_output}")

print(f"\n<<<{formatted_output}>>>")