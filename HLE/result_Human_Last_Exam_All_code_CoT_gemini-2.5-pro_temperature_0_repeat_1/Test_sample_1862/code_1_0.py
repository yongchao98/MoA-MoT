# Define the parameters for the initial scenario
d1_nm = 50  # Diameter of the first particle in nm
t1_s = 20   # Lithiation time for the first particle in seconds

# Define the parameters for the second scenario
d2_nm = 250 # Diameter of the second particle in nm

# Calculate the radius for each particle
# The characteristic diffusion length is the radius.
r1_nm = d1_nm / 2
r2_nm = d2_nm / 2

# The lithiation time (t) is proportional to the square of the radius (r).
# t_new / t_old = (r_new / r_old)^2
# t_new = t_old * (r_new / r_old)^2
t2_s = t1_s * (r2_nm / r1_nm)**2

# Print the final calculation and result
print(f"The lithiation time for the second particle is calculated as follows:")
print(f"t_new = t_old * (r_new / r_old)^2")
print(f"t_new = {t1_s} s * ({int(r2_nm)} nm / {int(r1_nm)} nm)^2")
print(f"t_new = {t1_s} s * ({int(r2_nm / r1_nm)})^2")
print(f"t_new = {t1_s} s * {int((r2_nm / r1_nm)**2)}")
print(f"t_new = {int(t2_s)} s")
<<<500>>>