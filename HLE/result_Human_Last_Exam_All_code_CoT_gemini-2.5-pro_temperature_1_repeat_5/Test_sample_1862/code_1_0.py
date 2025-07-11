# Define the parameters for the first case (full sphere)
d1 = 50  # diameter in nm
t1 = 20  # lithiation time in seconds

# Calculate the radius for the first case
r1 = d1 / 2

# Define the parameters for the second case (quarter sphere)
d2 = 250 # diameter in nm

# Calculate the radius for the second case.
# The characteristic diffusion length for a quarter sphere is its radius,
# representing the longest path from the surface to the innermost point.
r2 = d2 / 2

# Calculate the lithiation time for the second case using the diffusion time relationship:
# t2 = t1 * (r2 / r1)^2
t2 = t1 * (r2 / r1)**2

# Print the final equation with all the values
print("The lithiation time is proportional to the square of the particle radius (t ∝ r²).")
print("We can find the new time (t2) using the following relationship:")
print(f"t2 = t1 * (r2 / r1)²")
print(f"t2 = {t1} s * ({r2} nm / {r1} nm)²")
print(f"t2 = {t1} s * ({r2/r1})²")
print(f"t2 = {t1} s * {(r2/r1)**2}")
print(f"The calculated lithiation time is: {t2} seconds")