# Step 1: Define the given parameters for the two cases.
# Case 1: Sphere
t1 = 20.0  # Lithiation time in seconds
d1 = 50.0  # Diameter in nm

# Case 2: Quarter sphere
d2 = 250.0 # Diameter in nm

# Step 2: Calculate the radius for each particle.
# The characteristic diffusion length (L) is the radius (r) for both shapes.
r1 = d1 / 2.0
r2 = d2 / 2.0

# Step 3: Calculate the lithiation time for the second case (t2).
# The relationship is t2 / t1 = (r2 / r1)^2
# So, t2 = t1 * (r2 / r1)^2
t2 = t1 * (r2 / r1)**2

# Step 4: Print the explanation and the final calculated value.
print("The lithiation time is governed by the solid-state diffusion law, where time (t) is proportional to the square of the diffusion length (radius, r).")
print("The relationship is: t2 = t1 * (r2 / r1)^2")
print(f"\nGiven values:")
print(f"Initial time (t1): {t1} s")
print(f"Initial radius (r1): {d1} nm / 2 = {r1} nm")
print(f"New radius (r2): {d2} nm / 2 = {r2} nm")
print("\nCalculation:")
print(f"t2 = {t1} * ({r2} / {r1})^2")
print(f"t2 = {t1} * ({r2 / r1})^2")
print(f"t2 = {t1} * ({(r2 / r1)**2})")
print(f"Final Lithiation Time (t2): {t2} seconds")
print("<<<500.0>>>")