import math

# Step 1: Define the given constants and parameters.
# The charge of the particle is e.
# The mass of the particle is m (electron mass).
# The charge of the magnetic monopole is e_m = e / 16.
# The magnetic dipole moment is mu = 25 * mu_B, where mu_B = e*hbar / (2*m*c).
# The particle's speed is v = c / 100.
# The scattering angle is theta = pi / 30.

# Step 2: Write down the formulas for the differential cross sections.
# For monopole scattering: dsigma_mono/dOmega = (e * e_m / (2 * p * c * sin(theta/2)**2))**2
# where p = m*v is the momentum.
# For dipole scattering (averaged over orientation): dsigma_dip/dOmega = (e**2 * mu**2 / (3 * c**2 * hbar**2)) * cot(theta/2)**2

# Step 3: Form the ratio R = (dsigma_mono/dOmega) / (dsigma_dip/dOmega).
# R = [ (e * e_m / (2 * p * c * sin(theta/2)**2))**2 ] / [ (e**2 * mu**2 / (3 * c**2 * hbar**2)) * (cos(theta/2)/sin(theta/2))**2 ]
# R = (e**2 * e_m**2 / (4 * p**2 * c**2 * sin(theta/2)**4)) * (3 * c**2 * hbar**2 * sin(theta/2)**2) / (e**2 * mu**2 * cos(theta/2)**2)
# R = (3 * hbar**2 * e_m**2) / (4 * p**2 * mu**2 * sin(theta/2)**2 * cos(theta/2)**2)

# Step 4: Simplify the ratio using the identity sin(theta) = 2*sin(theta/2)*cos(theta/2).
# R = (3 * hbar**2 * e_m**2) / (p**2 * mu**2 * sin(theta)**2)

# Step 5: Substitute the given values for e_m, mu, and p.
# e_m = e / 16  => e_m**2 = e**2 / 256
# mu = 25 * (e*hbar / (2*m*c)) => mu**2 = 625 * e**2 * hbar**2 / (4 * m**2 * c**2)
# p = m*v = m * (c / 100) => p**2 = m**2 * c**2 / 10000

# R = (3 * hbar**2 * (e**2 / 256)) / ((m**2 * c**2 / 10000) * (625 * e**2 * hbar**2 / (4 * m**2 * c**2)) * sin(theta)**2)
# Many terms cancel out (e, hbar, m, c):
# R = (3 / 256) / ((1 / 10000) * (625 / 4) * sin(theta)**2)
# R = (3 / 256) * (40000 / 625) / sin(theta)**2
# R = (3 / 256) * 64 / sin(theta)**2
# R = (3 * 64 / 256) / sin(theta)**2
# R = 192 / 256 / sin(theta)**2
# R = 3 / (4 * sin(theta)**2)

# Step 6: Calculate the final numerical value.
theta_val = math.pi / 30
ratio = 3 / (4 * math.sin(theta_val)**2)

print("The formula for the ratio of the differential cross-sections simplifies to:")
print("R = 3 / (4 * sin(theta)**2)")
print(f"For theta = pi/30, the equation is:")
print(f"R = 3 / (4 * sin(pi/30)**2) = 3 / (4 * {math.sin(theta_val)**2:.6f})")
print("\nThe final numerical value for the ratio is:")
print(ratio)

# The final answer in the requested format
print(f"<<<{ratio:.4f}>>>")