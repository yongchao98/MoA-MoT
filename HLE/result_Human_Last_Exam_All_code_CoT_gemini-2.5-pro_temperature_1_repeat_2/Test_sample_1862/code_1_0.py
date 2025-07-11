import math

# --- Given parameters ---
# Case 1: Full sphere
diameter1_nm = 50
time1_s = 20

# Case 2: Quarter sphere
diameter2_nm = 250

# --- Calculations ---
# The characteristic diffusion length for both a full sphere and a quarter sphere
# is the radius of the original sphere.
radius1_nm = diameter1_nm / 2
radius2_nm = diameter2_nm / 2

# The lithiation time (t) is proportional to the square of the diffusion length (radius, r).
# t2 / t1 = (r2 / r1)^2
# Therefore, t2 = t1 * (r2 / r1)^2
time2_s = time1_s * (radius2_nm / radius1_nm)**2

# --- Output the result ---
print(f"The lithiation time for the new particle is calculated using the diffusion time relationship.")
print(f"The new lithiation time (t2) is given by the equation: t2 = t1 * (r2 / r1)²")
print(f"Substituting the given values:")
print(f"t2 = {time1_s} s * ({radius2_nm} nm / {radius1_nm} nm)²")
# Showing the intermediate step of the ratio squared
ratio_squared = (radius2_nm / radius1_nm)**2
print(f"t2 = {time1_s} s * {int(ratio_squared)}")
# Showing the final calculation
print(f"t2 = {int(time2_s)} s")
print("\nTherefore, the lithiation time for the 250 nm quarter-sphere particle is 500 seconds.")
