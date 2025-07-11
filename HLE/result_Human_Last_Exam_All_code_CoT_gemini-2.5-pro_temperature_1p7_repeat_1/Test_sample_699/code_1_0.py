import math

# --- Parameters and Assumptions ---
# Basquin exponent from the problem
b = -0.09
# Assumed number of cycles at the endurance limit (N_e).
# This is a standard engineering assumption for materials like steel to get a numerical answer.
N_e = 1e6
# Loading profile: cycle fractions (f_i) and their corresponding stress ratios (r_i = σ_i / σ_e)
cycle_fractions = [0.70, 0.20, 0.10]
stress_ratios = [1.0, 1.1, 1.2]

# --- Calculation ---
# According to the Palmgren-Miner rule, N_f = 1 / Σ(f_i / N_i).
# From the Basquin relationship, N_i = N_e * r_i^(1/b).
# Substituting this into the Miner's rule gives:
# N_f = N_e / Σ(f_i * r_i^(-1/b))

# 1. Calculate the exponent term used in the damage calculation
exponent = -1 / b

# 2. Calculate each term in the denominator of the N_f formula
# Each term is f_i * r_i^(exponent)
damage_terms = [f * (r ** exponent) for f, r in zip(cycle_fractions, stress_ratios)]

# 3. Sum the damage terms to get the denominator
damage_sum = sum(damage_terms)

# 4. Calculate the final fatigue life in cycles
N_f = N_e / damage_sum


# --- Output ---
print("--- Fatigue Life Estimation using Palmgren-Miner Rule ---")
print("The total fatigue life (N_f) is calculated based on linear damage accumulation.")
print("The governing equation is derived as: N_f = N_e / [Σ (f_i * (σ_i/σ_e)^(-1/b))]\n")
print("Key parameters and assumptions used:")
print(f"  Basquin exponent, b = {b}")
print(f"  Exponent term, -1/b = {exponent:.4f}")
print(f"  Cycles at endurance limit (assumed), N_e = {int(N_e):,}\n")
print("The final equation is built using these values:")

# Using f-strings to build and print the equation with numbers
equation_str = (
    f"N_f = {int(N_e):,} / [({cycle_fractions[0]:.2f} * {stress_ratios[0]:.1f}^{exponent:.3f}) "
    f"+ ({cycle_fractions[1]:.2f} * {stress_ratios[1]:.1f}^{exponent:.3f}) "
    f"+ ({cycle_fractions[2]:.2f} * {stress_ratios[2]:.1f}^{exponent:.3f})]"
)
print(equation_str)

# Print the equation with the intermediate terms calculated
values_str = (
    f"N_f = {int(N_e):,} / [{damage_terms[0]:.4f} "
    f"+ {damage_terms[1]:.4f} "
    f"+ {damage_terms[2]:.4f}]"
)
print(values_str)

# Print the equation with the summed denominator
denominator_str = f"N_f = {int(N_e):,} / {damage_sum:.4f}"
print(denominator_str)

# Print the final result
final_result_str = f"N_f = {N_f:,.0f} cycles\n"
print(final_result_str)

# Format the answer as requested: Aeb with 3 decimal places
formatted_answer = f"{N_f:.3e}".replace('e+0', 'e').replace('e+', 'e').replace('e-0', 'e-')
print(f"The estimated fatigue life in the required format is: {formatted_answer}")
print(f"<<<{formatted_answer}>>>")