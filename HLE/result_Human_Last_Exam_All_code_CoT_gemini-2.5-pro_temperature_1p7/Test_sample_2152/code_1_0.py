import math

# Given values
# Ratio of magnetic monopole charge to elementary charge
em_over_e = 1/16
# Dipole magnetic moment in units of Bohr magnetons
mu_over_muB = 25
# Particle speed as a fraction of the speed of light
v_over_c = 1/100
# Scattering angle in radians
theta = math.pi / 30

# The derived formula for the ratio R
# R = 3 / (16 * cos^2(theta/2))

# Perform the calculation
cos_theta_half = math.cos(theta / 2)
ratio = 3 / (16 * cos_theta_half**2)

# Print the final equation with numbers and the result
print(f"The ratio of the differential cross-sections is given by the formula: R = 3 / (16 * cos^2(theta/2))")
print(f"Substituting the angle theta = pi/30 rad:")
print(f"R = 3 / (16 * cos^2(({theta:.4f})/2))")
print(f"R = 3 / (16 * ({cos_theta_half:.4f})^2)")
print(f"R = 3 / (16 * {cos_theta_half**2:.4f})")
print(f"The final calculated ratio is: {ratio}")

# Output the answer in the specified format
# The final result is the numerical value of the ratio
# final_answer = ratio