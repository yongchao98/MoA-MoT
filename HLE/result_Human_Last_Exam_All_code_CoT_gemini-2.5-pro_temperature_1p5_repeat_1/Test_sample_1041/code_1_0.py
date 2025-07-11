import math

# Define constants
R = 1000.0
z0 = (0, 300)
w1 = (0, 0)
w2 = (2, 0)
gamma = 0.5772156649  # Euler-Mascheroni constant

# --- Numerator Calculation ---
# The numerator is g_D(z0, w1) + g_D(z0, w2)
# Distance from start z0 to w1
dist_z0_w1 = math.sqrt(z0[0]**2 + z0[1]**2)
# Distance from start z0 to w2
dist_z0_w2 = math.sqrt((z0[0] - w2[0])**2 + (z0[1] - w2[1])**2)

# g_D(z0, w1) approx (2/pi) * log(R / |z0 - w1|)
g_z0_w1 = (2 / math.pi) * math.log(R / dist_z0_w1)

# g_D(z0, w2) approx (2/pi) * log(R / |z0 - w2|)
g_z0_w2 = (2 / math.pi) * math.log(R / dist_z0_w2)

numerator = g_z0_w1 + g_z0_w2

# --- Denominator Calculation ---
# The denominator is g_D(w1, w1) + g_D(w1, w2)
# g_D(w1, w1) = g_D(0,0) approx (2/pi) * (log(R) + gamma + 2*log(2))
g_w1_w1 = (2 / math.pi) * (math.log(R) + gamma + 2 * math.log(2))

# g_D(w1, w2) = g_D(0,2), which is approx (2/pi) * log(R / |w1 - w2|)
dist_w1_w2 = math.sqrt((w1[0] - w2[0])**2 + (w1[1] - w2[1])**2)
g_w1_w2 = (2 / math.pi) * math.log(R / dist_w1_w2)

denominator = g_w1_w1 + g_w1_w2

# --- Final Probability ---
probability = numerator / denominator

# Print the calculation steps and the final result
print("--- Calculation of the hitting probability ---")
print(f"Starting point z0: {z0}")
print(f"Target set A: {{ {w1}, {w2} }}")
print(f"Disk radius R: {R}")
print("\n--- Intermediate Values (based on Green's function approximations) ---")
print("Numerator = g_D(z0, w1) + g_D(z0, w2)")
print(f"g_D(z0, w1) where |z0-w1|={dist_z0_w1:.4f}: {g_z0_w1:.4f}")
print(f"g_D(z0, w2) where |z0-w2|={dist_z0_w2:.4f}: {g_z0_w2:.4f}")
print(f"Total Numerator = {g_z0_w1:.4f} + {g_z0_w2:.4f} = {numerator:.4f}")

print("\nDenominator = g_D(w1, w1) + g_D(w1, w2)")
print(f"g_D(w1, w1) (visits to origin): {g_w1_w1:.4f}")
print(f"g_D(w1, w2) where |w1-w2|={dist_w1_w2:.4f}: {g_w1_w2:.4f}")
print(f"Total Denominator = {g_w1_w1:.4f} + {g_w1_w2:.4f} = {denominator:.4f}")

print("\n--- Final Result ---")
print(f"Probability = Numerator / Denominator")
print(f"Probability = {numerator:.4f} / {denominator:.4f} = {probability:.4f}")
print(f"\nThe probability rounded to three significant digits is {probability:.3g}")
<<<0.160>>>