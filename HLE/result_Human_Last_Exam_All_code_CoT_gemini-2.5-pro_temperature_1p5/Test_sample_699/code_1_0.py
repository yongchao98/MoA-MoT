import math

# --- 1. Define constants and assumptions ---
# Basquin exponent
b = -0.09
# Assumed number of cycles at the fatigue endurance limit (Ne)
N_e = 1_000_000
# Cycle fractions for each stress level
f1 = 0.70  # at 1.0 * sigma_e
f2 = 0.20  # at 1.1 * sigma_e
f3 = 0.10  # at 1.2 * sigma_e
# Stress ratios relative to the endurance limit (sigma_e)
s1 = 1.0
s2 = 1.1
s3 = 1.2

# --- 2. Calculate fatigue life (N_fi) for each stress level ---
# Using the formula N_fi = N_e * (sigma_i / sigma_e)^(1/b)
N_f1 = N_e * math.pow(s1, 1/b)
N_f2 = N_e * math.pow(s2, 1/b)
N_f3 = N_e * math.pow(s3, 1/b)

# --- 3. Apply Palmgren-Miner rule to find total life (N_T) ---
# N_T = 1 / (f1/N_f1 + f2/N_f2 + f3/N_f3)
# This simplifies to: N_T = N_e / (f1*(s1)^(-1/b) + f2*(s2)^(-1/b) + f3*(s3)^(-1/b))
exponent = -1/b
term1 = f1 * math.pow(s1, exponent)
term2 = f2 * math.pow(s2, exponent)
term3 = f3 * math.pow(s3, exponent)
denominator = term1 + term2 + term3
total_life_N_T = N_e / denominator

# --- 4. Print the calculation steps and final answer ---
print("Estimating fatigue life using the Palmgren-Miner linear damage rule.")
print("\nThe final life N_T is calculated as:")
print(f"N_T = N_e / (f1 * ({s1})^(-1/b) + f2 * ({s2})^(-1/b) + f3 * ({s3})^(-1/b))\n")

print("Given values and intermediate calculations:")
print(f"N_e = {int(N_e):,}")
print(f"b = {b}")
print(f"-1/b = {exponent:.4f}")

print("\nSubstituting the values into the equation:")
# We show the expanded denominator calculation
print(f"Denominator = ({f1:.2f} * {s1:.1f}^{exponent:.4f}) + ({f2:.2f} * {s2:.1f}^{exponent:.4f}) + ({f3:.2f} * {s3:.1f}^{exponent:.4f})")
print(f"Denominator = {term1:.4f} + {term2:.4f} + {term3:.4f}")
print(f"Denominator = {denominator:.4f}\n")

# Show the final equation with numbers
print("Final Equation:")
print(f"N_T = {int(N_e)} / {denominator:.4f}")

print(f"\nEstimated Fatigue Life (N_T) = {total_life_N_T:.3f} cycles")
final_answer = f"{total_life_N_T:.3e}"
print(f"Final answer in scientific notation (A.dddeB): {final_answer}")

print(f"\n<<<{final_answer}>>>")