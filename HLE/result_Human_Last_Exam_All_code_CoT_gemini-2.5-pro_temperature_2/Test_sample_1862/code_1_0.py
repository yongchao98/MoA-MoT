import math

# Initial particle parameters
d1_nm = 50  # diameter in nm
t1_s = 20   # time in seconds

# Second particle parameters
d2_nm = 250 # diameter in nm

# Calculate radii from diameters
r1_nm = d1_nm / 2
r2_nm = d2_nm / 2

# The lithiation time (t) is proportional to the square of the diffusion length (radius, r).
# t = k * r^2
# So, t2 / t1 = (r2^2) / (r1^2)
# t2 = t1 * (r2 / r1)^2

t2_s = t1_s * (r2_nm / r1_nm)**2

print("The calculation for the new lithiation time is based on the principle that diffusion time is proportional to the square of the diffusion length (the particle radius).")
print(f"Initial particle radius (r1) = {d1_nm} nm / 2 = {r1_nm} nm")
print(f"Initial lithiation time (t1) = {t1_s} s")
print(f"New particle radius (r2) = {d2_nm} nm / 2 = {r2_nm} nm")
print("\nFinal Equation:")
print(f"New Time = Initial Time * (New Radius / Initial Radius)^2")
print(f"New Time = {t1_s} s * ({r2_nm} nm / {r1_nm} nm)^2 = {t1_s} * ({r2_nm/r1_nm})^2 = {t2_s} seconds")
print(f"The calculated lithiation time for the quarter-sphere particle is {t2_s} seconds.")
<<<500>>>