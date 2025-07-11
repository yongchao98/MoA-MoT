import math

# --- Given Parameters ---
# Basquin exponent
b = -0.09
# Life fractions for each stress level
f1 = 0.70  # at sigma_e
f2 = 0.20  # at 1.1 * sigma_e
f3 = 0.10  # at 1.2 * sigma_e
# Stress level ratios relative to sigma_e
r2 = 1.1
r3 = 1.2

# --- Step 1: Calculate N1 (cycles to failure at endurance limit, sigma_e) ---
# This is based on the assumption that sigma_e = 0.5 * sigma_uts
# Derived from: sigma_e / sigma_uts = (N1 / 0.5)**b
# which leads to: N1 = 0.5 * (0.5)**(1/b)
N1 = 0.5 * (0.5)**(1/b)

# --- Step 2: Set up and solve the Palmgren-Miner equation for total life L ---
# The equation is: L = N1 / (f1 + f2*(r2**(-1/b)) + f3*(r3**(-1/b)))
# Calculate the denominator terms first for clarity.
term1 = f1
term2 = f2 * r2**(-1/b)
term3 = f3 * r3**(-1/b)
denominator = term1 + term2 + term3

# Calculate the final life L
L = N1 / denominator

# --- Step 3: Print the full calculation and the final answer ---
print("--- Fatigue Life Estimation ---")
print(f"Basquin exponent (b): {b}")
print(f"Assumed cycles to failure at endurance limit (N1): {N1:.3f}\n")

print("The total life (L) is calculated using the Palmgren-Miner rule:")
print(f"L = N1 / (f1 + f2*r2^(-1/b) + f3*r3^(-1/b))")
print("\nSubstituting the numerical values into the equation:")
print(f"L = {N1:.3f} / ({f1} + {f2} * {r2}^(-1/{b}) + {f3} * {r3}^(-1/{b}))")
print(f"L = {N1:.3f} / ({term1:.3f} + {term2:.3f} + {term3:.3f})")
print(f"L = {N1:.3f} / {denominator:.3f}")
print(f"L = {L:.3f} cycles\n")

print("Final Answer:")
print(f"Estimated Fatigue Life = {L:.3e} cycles")
