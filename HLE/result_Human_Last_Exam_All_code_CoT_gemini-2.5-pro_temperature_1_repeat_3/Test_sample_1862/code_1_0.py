import math

# Case 1: Initial small spherical particle
d1_nm = 50.0  # Diameter in nm
t1_s = 20.0   # Lithiation time in seconds

# Case 2: Final larger quarter-sphere particle
d2_nm = 250.0 # Diameter in nm

# The characteristic diffusion length (L) for both a full sphere and a quarter-sphere
# is the radius (r), as this is the longest path for Li-ions to travel.
r1_nm = d1_nm / 2.0
r2_nm = d2_nm / 2.0

# The lithiation time (t) is proportional to the square of the diffusion length (r^2).
# We can set up a ratio: t2 / t1 = (r2 / r1)^2
# We solve for t2: t2 = t1 * (r2 / r1)^2

t2_s = t1_s * (r2_nm / r1_nm)**2

print("The lithiation time is proportional to the square of the particle's radius.")
print("We can calculate the new time (t2) using the following relationship:")
print("t2 = t1 * (r2 / r1)^2\n")

print("Given values:")
print(f"Initial time (t1): {t1_s} s")
print(f"Initial radius (r1): {d1_nm} nm / 2 = {r1_nm} nm")
print(f"Final radius (r2): {d2_nm} nm / 2 = {r2_nm} nm\n")

print("Calculation:")
# The user prompt requested to output each number in the final equation.
# So we will print the equation with the numbers substituted.
print(f"t2 = {t1_s} s * ({r2_nm} nm / {r1_nm} nm)^2")
print(f"t2 = {t1_s} s * ({r2_nm/r1_nm})^2")
print(f"t2 = {t1_s} s * { (r2_nm/r1_nm)**2 }")
print(f"The new lithiation time is: {t2_s} seconds.")

<<<500.0>>>