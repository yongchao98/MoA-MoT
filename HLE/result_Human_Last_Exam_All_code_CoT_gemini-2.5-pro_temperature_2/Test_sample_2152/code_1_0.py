import math

# Plan:
# 1. Define the differential cross-section for a charged particle scattering from a magnetic monopole.
#    (d_sigma/d_Omega)_mono = (e * e_m / (m * v * c))^2 / sin(theta/2)^4
# 2. Define the differential cross-section for a charged particle scattering from a magnetic dipole
#    (in the first Born approximation, averaged over dipole orientations).
#    (d_sigma/d_Omega)_dip = (e^2 * mu^2) / (3 * h_bar^2 * c^2 * sin(theta/2)^2)
# 3. Compute the ratio R = (d_sigma/d_Omega)_mono / (d_sigma/d_Omega)_dip.
#    R = (3 * e_m^2 * h_bar^2) / (m^2 * v^2 * mu^2 * sin(theta/2)^2)
# 4. Substitute the given values:
#    e_m = e / 16
#    mu = 25 * mu_B = 25 * (e * h_bar) / (2 * m * c) (assuming the particle is an electron, m=m_e)
#    v = c / 100
# 5. Simplify the expression for R.
#    R = 3 / (4 * sin(theta/2)^2)
# 6. Calculate the numerical value.

# Scattering angle in radians
theta = math.pi / 30.0

# Calculate theta/2
theta_over_2 = theta / 2.0

# Calculate sin(theta/2)
sin_val = math.sin(theta_over_2)

# The final simplified formula for the ratio
# Ratio = 3 / (4 * sin^2(pi/60))
numerator = 3.0
denominator_coeff = 4.0
ratio = numerator / (denominator_coeff * sin_val**2)

print("The simplified expression for the ratio of the differential cross-sections (monopole/dipole) is:")
print("R = 3 / (4 * sin^2(theta/2))")
print("\nWith the given values:")
print(f"theta = pi/30 radians")
print(f"theta/2 = pi/60 radians = {theta_over_2:.5f} radians")
print("\nThe final equation with numerical values is:")
print(f"R = {numerator} / ({denominator_coeff} * sin({theta_over_2:.5f})^2)")
print(f"R = {numerator} / ({denominator_coeff} * {sin_val:.5f}^2)")
print(f"R = {numerator} / ({denominator_coeff} * {sin_val**2:.5f})")
print(f"R = {numerator} / ({denominator_coeff * sin_val**2:.5f})")
print(f"\nThe calculated ratio is: {ratio}")
