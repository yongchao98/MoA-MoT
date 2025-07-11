# Step 1: Define the initial parameters for the first particle (full sphere).
d1_nm = 50  # Diameter of the first particle in nm
t1_s = 20   # Lithiation time for the first particle in seconds

# Step 2: Define the parameters for the second particle (quarter sphere).
d2_nm = 250 # Diameter of the second particle in nm

# Step 3: Calculate the radius for each particle. The characteristic diffusion length is the radius.
r1_nm = d1_nm / 2
r2_nm = d2_nm / 2

# Step 4: The lithiation time (t) is proportional to the square of the diffusion length (radius, r).
# So, t2 / t1 = (r2 / r1)^2. We can calculate t2 from this relationship.
# t2 = t1 * (r2 / r1)^2

# Calculate the ratio of the radii
radius_ratio = r2_nm / r1_nm

# Calculate the new lithiation time
t2_s = t1_s * (radius_ratio ** 2)

# Step 5: Print the final calculation step-by-step.
print("The lithiation time is proportional to the square of the particle's radius (t ∝ r^2).")
print("The new lithiation time (t2) can be calculated using the formula: t2 = t1 * (r2 / r1)^2\n")
print(f"Given values:")
print(f"t1 = {t1_s} s")
print(f"r1 = {d1_nm} nm / 2 = {r1_nm} nm")
print(f"r2 = {d2_nm} nm / 2 = {r2_nm} nm\n")
print("Calculation:")
print(f"t2 = {t1_s} * ({r2_nm} / {r1_nm})^2")
print(f"t2 = {t1_s} * ({int(radius_ratio)})^2")
print(f"t2 = {t1_s} * {int(radius_ratio**2)}")
print(f"t2 = {int(t2_s)} seconds")

<<<500>>>