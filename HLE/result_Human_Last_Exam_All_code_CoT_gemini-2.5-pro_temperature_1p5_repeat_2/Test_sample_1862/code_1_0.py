import math

# Initial parameters for the first spherical particle
d1_nm = 50.0  # Diameter in nm
t1_s = 20.0   # Lithiation time in seconds

# Parameters for the second quarter-sphere particle
d2_nm = 250.0 # Diameter in nm

# The characteristic diffusion length for both shapes is the radius.
r1_nm = d1_nm / 2
r2_nm = d2_nm / 2

# The lithiation time (t) is proportional to the square of the radius (r).
# t2 / t1 = (r2^2) / (r1^2)
# t2 = t1 * (r2 / r1)^2
t2_s = t1_s * (r2_nm / r1_nm)**2

# Print the final equation with the calculated numbers
print(f"The new lithiation time (t2) is calculated using the formula: t2 = t1 * (r2 / r1)^2")
print(f"Substituting the values into the equation:")
print(f"t2 = {int(t1_s)} * ({int(r2_nm)} / {int(r1_nm)})^2")
print(f"t2 = {int(t1_s)} * ({int(r2_nm / r1_nm)})^2")
print(f"t2 = {int(t1_s)} * {int((r2_nm / r1_nm)**2)}")
print(f"Final calculated lithiation time: {int(t2_s)} seconds")

<<<500>>>