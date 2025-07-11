# Initial conditions
d1_nm = 50  # Diameter of the first particle in nm
t1_s = 20   # Lithiation time for the first particle in seconds

# Conditions for the second particle
d2_nm = 250 # Diameter of the second particle in nm

# Calculate radii
r1_nm = d1_nm / 2
r2_nm = d2_nm / 2

# The lithiation time (t) is proportional to the square of the radius (r).
# t2 = t1 * (r2 / r1)^2
t2_s = t1_s * (r2_nm / r1_nm)**2

# Output the equation and the final result
print(f"The initial particle has a radius of {r1_nm} nm and takes {t1_s} s to lithiate.")
print(f"The new particle has a radius of {r2_nm} nm.")
print("\nThe lithiation time scales with the square of the radius (diffusion length).")
print("The calculation is: t_new = t_initial * (r_new / r_initial)^2")
print("\nPlugging in the values:")
print(f"t_new = {t1_s} s * ({int(r2_nm)} nm / {int(r1_nm)} nm)^2")
print(f"t_new = {t1_s} s * ({int(r2_nm / r1_nm)})^2")
print(f"t_new = {t1_s} s * {int((r2_nm / r1_nm)**2)}")
print(f"Final lithiation time = {int(t2_s)} seconds")
