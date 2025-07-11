import math

# Define the velocity of particle A in the lab frame.
beta = 0.95

# Calculate the Lorentz factor, gamma.
gamma = 1 / math.sqrt(1 - beta**2)

# From the four-momentum of C in A's rest frame, we determine the emission angle theta*.
# C*_mu = (E_C*, P_C*/sqrt(2), 0, P_C*/sqrt(2))
# p*_z = P_C*/sqrt(2) and the perpendicular momentum p*_perp = P_C*/sqrt(2).
# The magnitude of the 3-momentum is |p*_C| = P_C*.
# cos(theta*) = p*_z / |p*_C| = 1/sqrt(2)
# sin(theta*) = p*_perp / |p*_C| = 1/sqrt(2)
cos_theta_star = 1 / math.sqrt(2)
sin_theta_star = 1 / math.sqrt(2)

# Use the Lorentz transformation formula for angles to find the angle theta in the lab frame.
# tan(theta) = sin(theta*) / (gamma * (cos(theta*) + beta))
tan_theta_lab = sin_theta_star / (gamma * (cos_theta_star + beta))

# Calculate the angle in degrees from its tangent.
theta_lab_rad = math.atan(tan_theta_lab)
theta_lab_deg = math.degrees(theta_lab_rad)

# Print the final equation with the calculated numbers and the result.
print("The angle θ in the lab frame is found using the equation: tan(θ) = sin(θ*)/(γ * (cos(θ*) + β))")
print("\nSubstituting the numerical values:")
print(f"tan(θ) = {sin_theta_star:.4f} / ({gamma:.4f} * ({cos_theta_star:.4f} + {beta}))")
print(f"tan(θ) = {sin_theta_star:.4f} / ({gamma:.4f} * {cos_theta_star + beta:.4f})")
print(f"tan(θ) = {sin_theta_star:.4f} / {gamma * (cos_theta_star + beta):.4f}")
print(f"tan(θ) = {tan_theta_lab:.4f}")
print(f"θ = arctan({tan_theta_lab:.4f})")

print(f"\nThe final angle is {theta_lab_deg:.3f} degrees.")
<<<7.586>>>