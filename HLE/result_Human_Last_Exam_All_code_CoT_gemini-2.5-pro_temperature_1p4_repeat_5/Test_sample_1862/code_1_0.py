# --- Input Parameters ---

# Particle 1: Full Sphere
d1 = 50.0  # nm, diameter
t1 = 20.0  # s, lithiation time

# Particle 2: Quarter Sphere
d2 = 250.0 # nm, diameter

# --- Calculations ---

# Step 1: Calculate radii
r1 = d1 / 2.0
r2 = d2 / 2.0

# Step 2: The lithiation time 't' is proportional to the characteristic length (Volume/Area).
# t_new / t_initial = (V2/A2) / (V1/A1)
# For a full sphere, V1/A1 = r1/3
# For a quarter sphere, V2/A2 = r2/6
# So, t2 = t1 * (r2/6) / (r1/3) which simplifies to t2 = t1 * (r2/r1) * 0.5

radius_ratio = r2 / r1
t2 = t1 * radius_ratio * 0.5

# --- Output ---

print("Calculating the new lithiation time (t2) for the quarter sphere.")
print("The calculation is based on the formula: t2 = t1 * (Radius_2 / Radius_1) * 0.5")
print("\n--- Final Equation ---")
print(f"t2 = {t1} s * ({r2} nm / {r1} nm) * 0.5")
print(f"t2 = {t1} s * {radius_ratio} * 0.5")
print(f"t2 = {t1 * radius_ratio} s * 0.5")
print(f"t2 = {t2} s")

# Final answer in the specified format is printed below this informative block.
# Do not copy this line. The final answer is at the very end.
# print(f"\n<<< {t2} >>>")

<<<50.0>>>