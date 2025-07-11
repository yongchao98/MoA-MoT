import math

# --- 1. Define constants and assumptions ---
b = -0.09  # Basquin exponent

# Loading cycle fractions
f_e = 0.70  # Fraction of life at endurance limit
f_1 = 0.20  # Fraction of life at 1.1 * sigma_e
f_2 = 0.10  # Fraction of life at 1.2 * sigma_e

# Stress ratios relative to the endurance limit
stress_ratio_1 = 1.1
stress_ratio_2 = 1.2

# Assumption: The fatigue endurance limit (sigma_e) is the stress
# that causes failure at N_e = 1,000,000 cycles.
N_e = 1e6

# --- 2. Calculate cycles to failure (N_i) for each stress level ---
# Using the formula: N_i = N_e * (stress_ratio_i)^(1/b)
N_1 = N_e * (stress_ratio_1)**(1/b)
N_2 = N_e * (stress_ratio_2)**(1/b)

# --- 3. Apply Palmgren-Miner rule to find total life (N_total) ---
# The rule is: (f_e * N_total / N_e) + (f_1 * N_total / N_1) + (f_2 * N_total / N_2) = 1
# Solving for N_total: N_total = 1 / (f_e/N_e + f_1/N_1 + f_2/N_2)
damage_sum_denominator = (f_e / N_e) + (f_1 / N_1) + (f_2 / N_2)
N_total = 1 / damage_sum_denominator

# --- 4. Print the results and the final equation ---
print("--- Fatigue Life Calculation ---")
print(f"Basquin Exponent (b): {b}")
print(f"Assumed life at endurance limit (N_e): {N_e:,.0f} cycles")
print("\nCalculating life at higher stress levels:")
print(f"Life at {stress_ratio_1}*sigma_e (N_1): {N_1:,.0f} cycles")
print(f"Life at {stress_ratio_2}*sigma_e (N_2): {N_2:,.0f} cycles")

print("\nUsing the Palmgren-Miner rule, the total life (N_total) is found by solving:")
print(f"1 = ({f_e} * N_total / {N_e:.0f}) + ({f_1} * N_total / {N_1:.0f}) + ({f_2} * N_total / {N_2:.0f})")
print("\nWhich simplifies to:")
print(f"N_total = 1 / ({f_e}/{N_e:.0f} + {f_1}/{N_1:.0f} + {f_2}/{N_2:.0f})")

print("\n--- Final Estimated Fatigue Life ---")
print(f"The total estimated fatigue life is {N_total:,.0f} cycles.")
print(f"In scientific notation (to 3 decimal places): {N_total:.3e}")

# The final answer in the requested format
final_answer = f"{N_total:.3e}"
# The prompt format is Aeb but python outputs Ae+b, so we remove the +
final_answer = final_answer.replace('e+', 'e')
print(f"\n<<<{final_answer}>>>")