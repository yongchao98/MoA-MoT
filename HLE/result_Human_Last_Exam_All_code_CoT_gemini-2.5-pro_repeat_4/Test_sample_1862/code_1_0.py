# Initial parameters from the problem description
d1_nm = 50.0  # Diameter of the first particle in nm
t1_s = 20.0   # Lithiation time for the first particle in seconds
d2_nm = 250.0 # Diameter of the second particle in nm

# Step 1: Calculate the radius for both particles
r1_nm = d1_nm / 2
r2_nm = d2_nm / 2

# Step 2: Calculate the ratio of the radii squared
radius_ratio_squared = (r2_nm / r1_nm)**2

# Step 3: Calculate the new lithiation time (t2)
# The lithiation time (t) is proportional to the radius (r) squared (t ∝ r²).
# Therefore, t2 = t1 * (r2/r1)²
t2_s = t1_s * radius_ratio_squared

# Step 4: Print the explanation and the step-by-step calculation
print("The lithiation time is proportional to the square of the particle's radius.")
print("The formula is: t_new = t_initial * (r_new / r_initial)²")
print("\n--- Calculation ---")
print(f"Initial time (t_initial): {t1_s} s")
print(f"Initial radius (r_initial): {d1_nm} nm / 2 = {r1_nm} nm")
print(f"New radius (r_new): {d2_nm} nm / 2 = {r2_nm} nm")
print("\nFinal Equation:")
print(f"t_new = {t1_s} s * ({r2_nm} nm / {r1_nm} nm)²")
print(f"t_new = {t1_s} s * ({r2_nm / r1_nm})²")
print(f"t_new = {t1_s} s * {radius_ratio_squared}")
print(f"t_new = {t2_s} s")

print(f"\nThe calculated lithiation time for the 250 nm quarter-sphere particle is {t2_s} seconds.")
print("<<<500.0>>>")