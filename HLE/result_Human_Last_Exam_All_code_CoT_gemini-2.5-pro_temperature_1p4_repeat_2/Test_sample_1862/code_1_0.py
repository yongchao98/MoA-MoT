# Initial parameters from the first case
t1 = 20  # seconds
d1 = 50  # nm

# Parameters for the second case
d2 = 250 # nm

# Step 1: Calculate the radius for each particle.
# The characteristic diffusion length is the radius.
r1 = d1 / 2
r2 = d2 / 2

# Step 2: Calculate the new lithiation time (t2).
# The lithiation time is proportional to the square of the radius (t ∝ r^2).
# Therefore, t2 = t1 * (r2/r1)^2
t2 = t1 * (r2 / r1)**2

# Step 3: Print the final calculation and the result.
# The formula is t2 = t1 * (r2 / r1)^2
print(f"The calculation for the new lithiation time (t2) is:")
print(f"t2 = {t1} s * ({r2} nm / {r1} nm)^2")
print(f"t2 = {t1} s * ({r2/r1})^2")
print(f"t2 = {t1} s * {(r2/r1)**2}")
print(f"t2 = {t2} s")

# Final Answer
# print(f"<<<{int(t2)}>>>")