import math

# Define the parameters for the initial spherical particle
d1_nm = 50.0  # Initial diameter in nm
t1_s = 20.0   # Initial lithiation time in seconds

# Define the parameters for the new quarter-sphere particle
d2_nm = 250.0 # New diameter in nm

# Calculate the radii from the diameters
r1_nm = d1_nm / 2.0
r2_nm = d2_nm / 2.0

# The lithiation process is limited by solid-state diffusion.
# The time required for diffusion (t) is proportional to the square of the characteristic diffusion length, which is the particle's radius (r).
# So, the relationship is: t_new / t_initial = (r_new / r_initial)^2
# We can solve for t_new: t_new = t_initial * (r_new / r_initial)^2

# Perform the calculation
t2_s = t1_s * (r2_nm / r1_nm)**2

# --- Output ---
print("The lithiation time (t) is proportional to the square of the particle radius (r), so t_new = t_initial * (r_new / r_initial)^2.")
print("\nCalculating the new lithiation time:")
print(f"Initial time (t_initial): {t1_s} s")
print(f"Initial radius (r_initial): {d1_nm} nm / 2 = {r1_nm} nm")
print(f"New radius (r_new): {d2_nm} nm / 2 = {r2_nm} nm")
print("\nFinal Equation:")
# The print statement below shows the full calculation with all the numbers.
print(f"t_new = {t1_s} s * ({r2_nm} nm / {r1_nm} nm)^2 = {t2_s:.0f} s")