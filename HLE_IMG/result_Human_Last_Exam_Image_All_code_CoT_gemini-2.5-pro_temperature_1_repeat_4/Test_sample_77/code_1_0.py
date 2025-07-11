import math

# This script calculates the particle size using the Debye-Scherrer equation.

# Step 1: Define the parameters from the problem description and the XRD graph.
# K is the Scherrer constant (dimensionless).
K = 0.9
# λ is the X-ray wavelength in nanometers.
lambda_nm = 0.1541
# 2θ is the position of the most intense peak (111) in degrees.
two_theta_deg = 38.7
# β is the estimated Full Width at Half Maximum (FWHM) of the peak in degrees.
# This is estimated visually from the graph for the peak at 38.7°.
fwhm_deg = 0.6

# Step 2: Convert angles from degrees to radians for the calculation.
# The Bragg angle (θ) is half of the 2θ value.
theta_rad = math.radians(two_theta_deg / 2)
# The FWHM (β) must also be in radians.
beta_rad = math.radians(fwhm_deg)

# Step 3: Calculate the particle size (D) using the Debye-Scherrer equation.
# D = (K * λ) / (β * cos(θ))
particle_size_nm = (K * lambda_nm) / (beta_rad * math.cos(theta_rad))

# Step 4: Print the parameters, the final equation with values, and the result.
print("--- Particle Size Estimation using Debye-Scherrer Equation ---")
print("Equation: D = (K * λ) / (β * cos(θ))")
print("\nParameters used for calculation:")
print(f"Scherrer Constant (K): {K}")
print(f"X-ray Wavelength (λ): {lambda_nm} nm")
print(f"Peak Position (2θ): {two_theta_deg}°")
print(f"Estimated FWHM (β): {fwhm_deg}°")

print("\nCalculation steps:")
print(f"1. Convert Bragg angle to radians: θ = ({two_theta_deg} / 2)° = {math.degrees(theta_rad):.2f}° = {theta_rad:.4f} rad")
print(f"2. Convert FWHM to radians: β = {fwhm_deg}° = {beta_rad:.4f} rad")
print(f"3. Calculate cos(θ): cos({theta_rad:.4f}) = {math.cos(theta_rad):.4f}")

print("\nFinal Equation with substituted values:")
print(f"D = ({K} * {lambda_nm}) / ({beta_rad:.4f} * {math.cos(theta_rad):.4f})")
print(f"D = {K * lambda_nm:.4f} / {(beta_rad * math.cos(theta_rad)):.4f}")

print("\nResult:")
print(f"The estimated particle size (D) is {particle_size_nm:.2f} nm.")