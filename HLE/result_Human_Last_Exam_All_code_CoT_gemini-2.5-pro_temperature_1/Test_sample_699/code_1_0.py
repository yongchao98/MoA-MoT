import math

# --- Given parameters and assumptions ---
# Basquin exponent
b = -0.09
# The problem statement is missing a relationship between the endurance limit (sigma_e)
# and the ultimate tensile strength (sigma_uts). A common engineering assumption is made here.
sigma_e_over_sigma_uts = 0.5

# --- Step 1: Relate the S-N curve parameters ---
# From Basquin's law: sigma_a = C * N^b
# We have two points on the curve:
# 1) sigma_uts = C * (0.5)^b
# 2) sigma_e = C * (N_e)^b
# Dividing (2) by (1): sigma_e / sigma_uts = (N_e / 0.5)^b
# Solving for N_e (life at the endurance limit): N_e = 0.5 * (sigma_e / sigma_uts)^(1/b)

# --- Step 2: Calculate the life at each stress level ---
inv_b = 1.0 / b

# Calculate N1, the life at the endurance limit sigma_e
N1 = 0.5 * (sigma_e_over_sigma_uts)**inv_b

# Calculate N2, the life at stress level 1.1 * sigma_e
# From sigma_2 / sigma_1 = (N2 / N1)^b  =>  1.1 = (N2/N1)^b
# N2 = N1 * (1.1)^(1/b)
N2 = N1 * (1.1)**inv_b

# Calculate N3, the life at stress level 1.2 * sigma_e
# N3 = N1 * (1.2)^(1/b)
N3 = N1 * (1.2)**inv_b

# --- Step 3: Apply the Palmgren-Miner Rule ---
# The rule states: Sum(n_i / N_i) = 1
# (0.7 * N_f / N1) + (0.2 * N_f / N2) + (0.1 * N_f / N3) = 1
# Rearranging for N_f: 1 / N_f = (0.7 / N1) + (0.2 / N2) + (0.1 / N3)

# Calculate each term in the sum (damage per cycle for a load block of 1)
damage_fraction_1 = 0.7 / N1
damage_fraction_2 = 0.2 / N2
damage_fraction_3 = 0.1 / N3

# Total damage is the sum of the fractions
total_damage_per_cycle_block = damage_fraction_1 + damage_fraction_2 + damage_fraction_3

# The total number of cycles to failure (N_f) is the reciprocal of the total damage
N_f = 1.0 / total_damage_per_cycle_block

# --- Step 4: Display the results and the final equation ---
print(f"Solving for fatigue life using the Palmgren-Miner rule.")
print(f"An assumption is required: sigma_e / sigma_uts = {sigma_e_over_sigma_uts}\n")
print(f"Basquin exponent b = {b}")
print("--- Intermediate Life Calculations ---")
print(f"Life at sigma_e (N1) = {N1:.3f} cycles")
print(f"Life at 1.1*sigma_e (N2) = {N2:.3f} cycles")
print(f"Life at 1.2*sigma_e (N3) = {N3:.3f} cycles\n")

print("--- Palmgren-Miner Equation ---")
print(f"1/N_f = (Cycle fraction at level 1 / N1) + (Cycle fraction at level 2 / N2) + (Cycle fraction at level 3 / N3)")
print(f"1/N_f = (0.7 / {N1:.3f}) + (0.2 / {N2:.3f}) + (0.1 / {N3:.3f})")
print(f"1/N_f = {damage_fraction_1:.6f} + {damage_fraction_2:.6f} + {damage_fraction_3:.6f}")
print(f"1/N_f = {total_damage_per_cycle_block:.6f}")
print(f"N_f = 1 / {total_damage_per_cycle_block:.6f} = {N_f:.3f} cycles\n")

# --- Step 5: Format the final answer ---
if N_f != 0:
    exponent = math.floor(math.log10(abs(N_f)))
    mantissa = N_f / 10**exponent
    print("--- Final Answer ---")
    print(f"The estimated fatigue life is {mantissa:.3f}e{exponent} cycles.")
    final_answer_str = f"<<<{mantissa:.3f}e{exponent}>>>"
else:
    final_answer_str = "<<<0.000e0>>>"

print(final_answer_str)