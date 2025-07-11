# Case 1: Initial full sphere particle
d1_nm = 50  # diameter in nm
t1_s = 20   # lithiation time in seconds

# Case 2: Quarter sphere particle
d2_nm = 250 # diameter in nm

# The characteristic diffusion length (L) for a spherical geometry is its radius (r).
# This holds true even for a quarter sphere, as the longest diffusion path is still from the surface to the center point (the corner of the quarter sphere).
r1_nm = d1_nm / 2
r2_nm = d2_nm / 2

# The lithiation time (t) is proportional to the square of the diffusion length (r^2).
# t2 / t1 = (r2 / r1)^2
# Therefore, t2 = t1 * (r2 / r1)^2
t2_s = t1_s * (r2_nm / r1_nm)**2

# Print the final equation with all the numbers
print("The lithiation time for the second particle (t2) can be calculated using the relationship:")
print("t2 = t1 * (r2 / r1)^2")
print(f"t2 = {t1_s} s * ({r2_nm} nm / {r1_nm} nm)^2")
print(f"t2 = {t1_s} s * ({r2_nm/r1_nm})^2")
print(f"t2 = {t1_s} s * ({(r2_nm/r1_nm)**2})")
print(f"Calculated Lithiation Time: {t2_s} seconds")