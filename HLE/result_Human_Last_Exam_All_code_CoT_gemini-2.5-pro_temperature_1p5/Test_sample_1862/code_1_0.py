import math

# Case 1: Initial particle information
d1_nm = 50.0  # Diameter of the first particle in nm
t1_s = 20.0   # Lithiation time for the first particle in seconds

# Case 2: New particle information
d2_nm = 250.0 # Diameter of the second particle in nm

# The characteristic diffusion length (L) for a sphere or part of a sphere is its radius (r).
# The time for diffusion (t) is proportional to the square of the diffusion length (t ∝ L²).
# So, t2 / t1 = (L2 / L1)² = (r2 / r1)²
# t2 = t1 * (r2 / r1)²

# Calculate the radii
r1_nm = d1_nm / 2
r2_nm = d2_nm / 2

# Calculate the new lithiation time (t2)
t2_s = t1_s * (r2_nm / r1_nm)**2

# Print the final equation with all numbers included
print("The lithiation time (t2) is calculated using the diffusion time relationship:")
print(f"t2 = t1 * (r2 / r1)²")
print(f"t2 = {t1_s} s * ({r2_nm} nm / {r1_nm} nm)²")
print(f"t2 = {t1_s} s * ({r2_nm / r1_nm})²")
print(f"t2 = {t1_s} s * ({(r2_nm / r1_nm)**2})")
print(f"Final calculated lithiation time (t2): {t2_s} seconds")
<<<500.0>>>