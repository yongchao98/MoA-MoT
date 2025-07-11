import math

# --- Parameters for Case 1: Full Sphere ---
d1_nm = 50.0  # Initial diameter in nm
t1_s = 20.0   # Initial lithiation time in seconds

# --- Parameters for Case 2: Quarter Sphere ---
d2_nm = 250.0 # New diameter in nm

# --- Calculations ---

# Radii of the particles
r1_nm = d1_nm / 2
r2_nm = d2_nm / 2

# The lithiation time (t) is proportional to the Volume (V) and
# inversely proportional to the Surface Area (A).
# t ∝ V/A

# For the full sphere (Case 1), the ratio V1/A1 simplifies to r1/3.
# For the quarter sphere (Case 2), the ratio V2/A2 simplifies to r2/6.

# The relationship is t2 / t1 = (V2/A2) / (V1/A1)
# t2 = t1 * ( (r2/6) / (r1/3) )
# t2 = t1 * (r2/r1) * (3/6)
# t2 = t1 * (r2/r1) * 0.5

# Calculate the new lithiation time
t2_s = t1_s * (r2_nm / r1_nm) * 0.5

# --- Output ---
print("The lithiation time for the quarter-sphere particle is calculated based on the ratio of its volume to surface area.")
print("The formula derived from the physical model is: t_new = t_initial * (r_new / r_initial) * 0.5")
print("\nFinal Calculation:")
print(f"{t2_s:.1f} s = {t1_s:.1f} s * ({r2_nm:.1f} nm / {r1_nm:.1f} nm) * 0.5")
<<<50.0>>>