import math

# --- Given and Estimated Parameters ---
# Scherrer constant (dimensionless)
K = 0.9
# X-ray wavelength in nanometers (nm)
lambda_nm = 0.1541
# Peak position from XRD pattern in degrees (2*theta)
two_theta_deg = 38.7
# Full Width at Half Maximum (FWHM) of the peak in degrees, estimated from the graph
fwhm_deg = 0.6

# --- Calculations ---
# 1. Calculate the Bragg angle (theta) in radians
theta_rad = math.radians(two_theta_deg / 2)

# 2. Convert FWHM (beta) from degrees to radians
fwhm_rad = math.radians(fwhm_deg)

# 3. Calculate the crystallite size (D) using the Scherrer equation
# D = (K * lambda) / (beta * cos(theta))
particle_size_nm = (K * lambda_nm) / (fwhm_rad * math.cos(theta_rad))

# --- Output the results ---
print("--- Particle Size Estimation using Scherrer Equation ---")
print(f"Scherrer Equation: D = (K * \u03BB) / (\u03B2 * cos(\u03B8))")
print("\n--- Input Values ---")
print(f"Scherrer Constant (K): {K}")
print(f"X-ray Wavelength (\u03BB): {lambda_nm} nm")
print(f"Peak Position (2\u03B8): {two_theta_deg} degrees")
print(f"Bragg Angle (\u03B8): {two_theta_deg / 2:.2f} degrees")
print(f"Full Width at Half Maximum (\u03B2): {fwhm_deg} degrees")

print("\n--- Calculation ---")
print(f"D = ({K} * {lambda_nm}) / ({fwhm_rad:.5f} rad * cos({theta_rad:.5f} rad))")
print(f"D = {K * lambda_nm:.5f} / ({fwhm_rad:.5f} * {math.cos(theta_rad):.5f})")
print(f"D = {K * lambda_nm:.5f} / {(fwhm_rad * math.cos(theta_rad)):.5f}")

print("\n--- Result ---")
print(f"Estimated Particle Size (D): {particle_size_nm:.2f} nm")
