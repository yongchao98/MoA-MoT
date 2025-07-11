# Initial parameters for the first particle
d1_nm = 50  # Diameter in nm
t1_s = 20   # Time in seconds

# Parameters for the second particle
d2_nm = 250 # Diameter in nm

# The shape being a quarter-sphere does not change the characteristic diffusion length,
# which is the radius. The time for full lithiation is determined by the time it
# takes for lithium ions to diffuse to the furthest point from the surface (the center/corner),
# which is the radius.

# Calculate the radii
r1_nm = d1_nm / 2
r2_nm = d2_nm / 2

# The time for a diffusion-limited process is proportional to the square of the
# characteristic diffusion length (the radius in this case).
# t2 / t1 = (r2 / r1)^2
# Therefore, t2 = t1 * (r2 / r1)^2
t2_s = t1_s * (r2_nm / r1_nm)**2

# Print the final equation with all the numbers and the result
print(f"The lithiation time for the second particle is calculated as:")
print(f"Time = Initial Time * (Final Radius / Initial Radius)^2")
print(f"Time = {t1_s} s * ({int(r2_nm)} nm / {int(r1_nm)} nm)^2")
print(f"Time = {t1_s} s * ({int(r2_nm/r1_nm)})^2")
print(f"Time = {t1_s} s * {int((r2_nm/r1_nm)**2)}")
print(f"Time = {int(t2_s)} s")
