import math

# --- Step 1: Define the given parameters ---
# Basquin exponent
b = -0.09
# Life fractions for each stress level
f1, f2, f3 = 0.70, 0.20, 0.10
# Stress multipliers for each level (relative to the endurance limit)
m1, m2, m3 = 1.0, 1.1, 1.2

# --- Step 2: State the assumption for N1 ---
# The life at the endurance limit (N1) is not given.
# A standard engineering assumption for the endurance limit is 10^7 cycles.
N1 = 1e7

print("Step-by-Step Calculation:")
print(f"1. Assuming life at endurance limit (N1): {N1:.0e} cycles")

# --- Step 3: Calculate the Basquin life exponent for the denominator ---
# The exponent in the formula is -1/b
life_exp = -1.0 / b
print(f"2. Calculating the exponent term (-1/b): -1 / {b} = {life_exp:.4f}")

# --- Step 4: Calculate the terms in the denominator of the life equation ---
term1 = f1 * (m1 ** life_exp)  # m1 is 1.0, so this is just f1
term2 = f2 * (m2 ** life_exp)
term3 = f3 * (m3 ** life_exp)

print(f"3. Calculating the denominator of the Palmgren-Miner equation:")
print(f"   Denominator = ({f1} * {m1}**({life_exp:.2f})) + ({f2} * {m2}**({life_exp:.2f})) + ({f3} * {m3}**({life_exp:.2f}))")
print(f"   Term 1 (at 1.0*σ_e): {term1:.4f}")
print(f"   Term 2 (at 1.1*σ_e): {term2:.4f}")
print(f"   Term 3 (at 1.2*σ_e): {term3:.4f}")

# --- Step 5: Calculate the final life L ---
denominator = term1 + term2 + term3
total_life = N1 / denominator

print(f"4. Summing the denominator terms: {denominator:.4f}")
print(f"5. Calculating total life (L = N1 / Denominator):")
print(f"   L = {N1:.0e} / {denominator:.4f}")
print(f"   L = {total_life:.0f} cycles")

# --- Step 6: Format the final answer ---
# Format to 'Aeb' with 3 decimal places (e.g., 1.234e+06)
formatted_life = f"{total_life:.3e}"

print("\n--- Final Answer ---")
print("The estimated fatigue life of the specimen is:")
print(formatted_life)
<<<4.677e+06>>>