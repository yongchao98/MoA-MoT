import math

# Given parameters for the first case (full sphere)
d1_nm = 50.0  # diameter in nm
t1_s = 20.0   # lithiation time in seconds

# Given parameters for the second case (quarter sphere)
d2_nm = 250.0 # diameter in nm

# Step 1: Calculate the radii from the diameters
r1_nm = d1_nm / 2
r2_nm = d2_nm / 2

# Step 2: Apply the diffusion principle.
# The lithiation process is limited by solid-state diffusion. The characteristic time (t) for diffusion
# is proportional to the square of the characteristic diffusion length (L), which is the radius (r) for a sphere or a piece of a sphere.
# So, t ∝ r².
# We can set up a ratio: t2 / t1 = (r2 / r1)²

# Step 3: Solve for t2
# t2 = t1 * (r2 / r1)²
ratio_r = r2_nm / r1_nm
t2_s = t1_s * (ratio_r)**2

# Step 4: Print the explanation and the final calculation
print("The lithiation time (t) is proportional to the square of the particle's radius (r), based on diffusion principles (t ∝ r²).")
print(f"Given parameters:")
print(f"  - Time for the first particle (t1): {t1_s} s")
print(f"  - Radius of the first particle (r1): {d1_nm} nm / 2 = {r1_nm} nm")
print(f"  - Radius of the second particle (r2): {d2_nm} nm / 2 = {r2_nm} nm")
print("\nWe can calculate the new time (t2) using the following relationship:")
print("  t2 = t1 * (r2 / r1)²")
print("\nCalculation:")
# Displaying the calculation with all the numbers as requested.
print(f"  t2 = {t1_s} * ({r2_nm} / {r1_nm})²")
print(f"  t2 = {t1_s} * ({int(ratio_r)})²")
print(f"  t2 = {t1_s} * {int(ratio_r**2)}")
print(f"  t2 = {int(t2_s)} s")

print(f"\nThe calculated lithiation time for the 250 nm quarter-sphere particle is {int(t2_s)} seconds.")
print("<<<500>>>")