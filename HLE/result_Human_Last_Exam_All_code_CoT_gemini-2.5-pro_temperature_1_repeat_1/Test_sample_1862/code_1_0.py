import math

# Initial parameters for the first particle (full sphere)
d1_nm = 50  # diameter in nm
t1_s = 20   # lithiation time in seconds

# Parameters for the second particle (quarter sphere)
d2_nm = 250 # diameter in nm

# The rate-limiting step is solid-state diffusion.
# The lithiation time (t) is proportional to the square of the diffusion length,
# which is the radius (r) of the particle: t ∝ r^2.
# So, t2 / t1 = (r2 / r1)^2

# Calculate the radii
r1_nm = d1_nm / 2
r2_nm = d2_nm / 2

# Calculate the new lithiation time (t2)
# t2 = t1 * (r2 / r1)^2
t2_s = t1_s * (r2_nm / r1_nm)**2

print("Step 1: Define the initial conditions.")
print(f"  - First particle (sphere) diameter (d1): {d1_nm} nm")
print(f"  - First particle radius (r1): {d1_nm} / 2 = {r1_nm} nm")
print(f"  - First particle lithiation time (t1): {t1_s} s")
print(f"  - Second particle (quarter sphere) diameter (d2): {d2_nm} nm")
print(f"  - Second particle radius (r2): {d2_nm} / 2 = {r2_nm} nm\n")

print("Step 2: Apply the diffusion time relationship.")
print("The formula is: t2 = t1 * (r2 / r1)^2\n")

print("Step 3: Calculate the final lithiation time (t2).")
print(f"t2 = {t1_s} s * ({r2_nm} nm / {r1_nm} nm)^2")
print(f"t2 = {t1_s} s * ({int(r2_nm/r1_nm)})^2")
print(f"t2 = {t1_s} s * {int((r2_nm/r1_nm)**2)}")
print(f"t2 = {int(t2_s)} s\n")

print(f"The calculated lithiation time for the 250 nm quarter-sphere particle is {int(t2_s)} seconds.")
<<<500>>>