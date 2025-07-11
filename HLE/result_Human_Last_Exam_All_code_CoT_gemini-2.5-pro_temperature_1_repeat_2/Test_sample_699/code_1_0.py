import math

# Step 1: Define the given constants from the problem statement.
b = -0.09  # Basquin exponent

# Step 2: State the necessary assumption for the endurance ratio.
# This ratio is not provided, so a common engineering approximation is used.
# This is the ratio of the fatigue endurance limit to the ultimate tensile strength.
sigma_e_over_sigma_uts = 0.5

# Step 3: Calculate the number of cycles to failure at the endurance limit (N_e).
# From Basquin's law, we have two points:
# (1) sigma_uts = C * (0.5)^b
# (2) sigma_e = C * (N_e)^b
# Dividing (2) by (1) gives: (sigma_e / sigma_uts) = (N_e / 0.5)^b
# Rearranging for N_e gives: N_e = 0.5 * (sigma_e / sigma_uts)^(1/b)
N_e = 0.5 * (sigma_e_over_sigma_uts)**(1/b)

# Step 4: Define the variables for the Palmgren-Miner rule.
# The total life N_life is subjected to a cycle distribution:
# n1 = 0.7 * N_life at sigma_1 = sigma_e
# n2 = 0.2 * N_life at sigma_2 = 1.1 * sigma_e
# n3 = 0.1 * N_life at sigma_3 = 1.2 * sigma_e
# Miner's Rule: (n1/N_f1) + (n2/N_f2) + (n3/N_f3) = 1
# Substituting N_fi = N_e * (sigma_e/sigma_i)^(-1/b) and n_i leads to:
# N_life = N_e / (0.7 + 0.2 * 1.1^(-1/b) + 0.1 * 1.2^(-1/b))

# Step 5: Calculate the terms in the denominator of the life equation.
inv_b_pos = -1/b
term2_multiplier = 1.1**inv_b_pos
term3_multiplier = 1.2**inv_b_pos
denominator = 0.7 + 0.2 * term2_multiplier + 0.1 * term3_multiplier

# Step 6: Calculate the total estimated fatigue life.
N_life = N_e / denominator

# Step 7: Output the numbers used in the final equation and the result.
print("Calculation of the total fatigue life (N_life):")
print(f"The equation for total life is: N_life = N_e / (0.7 + 0.2 * 1.1^(-1/b) + 0.1 * 1.2^(-1/b))")
print("\n--- Intermediate Values ---")
print(f"Basquin exponent (b): {b}")
print(f"Life at endurance limit (N_e), assuming sigma_e/sigma_uts = {sigma_e_over_sigma_uts}: {N_e:.3f} cycles")
print(f"Value of -1/b: {inv_b_pos:.3f}")
print(f"Term for 1.1*sigma_e stress level (1.1^(-1/b)): {term2_multiplier:.3f}")
print(f"Term for 1.2*sigma_e stress level (1.2^(-1/b)): {term3_multiplier:.3f}")
print(f"Final denominator value: {denominator:.3f}")

print("\n--- Final Equation with Numbers ---")
# The following line shows the complete equation with all the calculated numerical values plugged in.
print(f"N_life = {N_e:.3f} / (0.7 + 0.2 * {term2_multiplier:.3f} + 0.1 * {term3_multiplier:.3f})")
print(f"N_life = {N_e:.3f} / {denominator:.3f}")

print("\n--- Estimated Fatigue Life ---")
# The final answer is formatted to 3 decimal places in scientific notation.
print(f"The estimated fatigue life is {N_life:.3e} cycles.")