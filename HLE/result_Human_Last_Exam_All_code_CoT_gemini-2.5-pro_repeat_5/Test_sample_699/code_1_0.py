import math

# --- Given parameters ---
b = -0.09  # Basquin exponent

# --- Step 1 & 2: Assumptions and defining life at endurance limit ---
# We assume the endurance limit sigma_e corresponds to a life of 10^7 cycles.
N1 = 1e7  # Life (cycles) at the endurance limit, sigma_e

# The problem states the load is composed of:
# 70% of life at sigma_e (Stress level 1)
# 20% of life at 1.1*sigma_e (Stress level 2)
# 10% of life at 1.2*sigma_e (Stress level 3)
# These are the cycle proportions in a repeating block of loading.
c1_prop = 0.70
c2_prop = 0.20
c3_prop = 0.10

# --- Step 3: Calculate life at each stress level ---
# From the Basquin relation, N = (sigma_a / C)^(1/b).
# If we have two stress levels sigma_a and sigma_b, the ratio of their lives is:
# N_a / N_b = (sigma_a / C)^(1/b) / (sigma_b / C)^(1/b) = (sigma_a / sigma_b)^(1/b)
# So, N_a = N_b * (sigma_a / sigma_b)^(1/b)

# We know N1 at sigma1 = sigma_e.
# For N2, sigma2 = 1.1 * sigma_e, so sigma2/sigma1 = 1.1
# For N3, sigma3 = 1.2 * sigma_e, so sigma3/sigma1 = 1.2

one_over_b = 1.0 / b
N2 = N1 * math.pow(1.1, one_over_b)
N3 = N1 * math.pow(1.2, one_over_b)

# --- Step 4 & 5: Apply Palmgren-Miner rule ---
# The Palmgren-Miner rule for total life N_f is:
# N_f = 1 / ( (c1_prop / N1) + (c2_prop / N2) + (c3_prop / N3) )
# Note: This is a simplified form assuming a block of 1 cycle, where fractions are applied.
# The result is the same as considering a block of 100 cycles (70, 20, 10).

# Calculate the total damage per cycle block (the denominator)
damage_per_block = (c1_prop / N1) + (c2_prop / N2) + (c3_prop / N3)

# The total life is the inverse of the damage per block
N_f = 1.0 / damage_per_block

# --- Step 6: Format and print the output ---
print("This script estimates fatigue life using the Palmgren-Miner rule.")
print("-" * 60)
print("The final calculation is based on the formula:")
print("N_f = 1 / [ (proportion_1 / N1) + (proportion_2 / N2) + (proportion_3 / N3) ]\n")
print("Where:")
print(f"Cycle Proportion 1 = {c1_prop}, Life N1 (at \u03C3_e)      = {N1:.1f} cycles")
print(f"Cycle Proportion 2 = {c2_prop}, Life N2 (at 1.1*\u03C3_e) = {N2:.1f} cycles")
print(f"Cycle Proportion 3 = {c3_prop}, Life N3 (at 1.2*\u03C3_e) = {N3:.1f} cycles")
print("-" * 60)

# Print the final equation with all numbers filled in
print("Final equation with values:")
print(f"N_f = 1 / [ ({c1_prop} / {N1:.1f}) + ({c2_prop} / {N2:.1f}) + ({c3_prop} / {N3:.1f}) ]\n")

# Print the final answer in the requested format
print(f"Estimated total fatigue life (N_f): {'{:.3e}'.format(N_f)} cycles")

# Final answer in specified format
final_answer_str = "{:.3e}".format(N_f)
print(f"\n<<< {final_answer_str} >>>")