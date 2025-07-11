# Step 1: Define the known variables from the problem statement.
d1_nm = 50  # Diameter of the first spherical particle in nm
t1_s = 20   # Lithiation time for the first particle in seconds
d2_nm = 250 # Diameter of the second quarter-sphere particle in nm

# Step 2: Calculate the radius for each particle. The characteristic diffusion length is the radius.
r1_nm = d1_nm / 2
r2_nm = d2_nm / 2

# Step 3: State the principle. The time (t) is proportional to the square of the radius (r).
# t = k * r^2, where k is a constant.
# Therefore, t2 / t1 = (r2^2) / (r1^2).
# We can solve for t2: t2 = t1 * (r2 / r1)^2

# Step 4: Calculate the ratio of the radii squared.
radius_ratio_squared = (r2_nm / r1_nm)**2

# Step 5: Calculate the new lithiation time (t2).
t2_s = t1_s * radius_ratio_squared

# Step 6: Print the final calculation and the result.
print("The lithiation time (t) is proportional to the square of the particle radius (r).")
print("t_new / t_initial = (r_new / r_initial)^2")
print(f"t_new = t_initial * (r_new / r_initial)^2")
print(f"t_new = {t1_s} s * ({r2_nm} nm / {r1_nm} nm)^2")
print(f"t_new = {t1_s} s * ({r2_nm/r1_nm})^2")
print(f"t_new = {t1_s} s * {radius_ratio_squared}")
print(f"The calculated new lithiation time is: {t2_s} seconds.")
<<<500.0>>>