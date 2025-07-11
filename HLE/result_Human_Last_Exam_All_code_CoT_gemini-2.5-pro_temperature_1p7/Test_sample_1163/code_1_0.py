import math

# Step 1: Define constants
period = 26000  # years
t_A = -3000     # years
t_B = 10000     # years
epsilon_deg = 23.5
epsilon_rad = math.radians(epsilon_deg)

# Step 2: Key insight from the problem statement
# The time between star A's last crossing and star B's next crossing is 13000 years,
# which is exactly half the precession period.
# This implies two geometric conditions on their ecliptic coordinates (lambda, beta):
# 1. Their ecliptic longitudes today are the same: lambda_A = lambda_B = lambda_0
# 2. Their ecliptic latitudes are opposite: beta_A = -beta_B = beta
# (A detailed derivation shows this is the only configuration consistent with
# both stars having the same crossing time pattern).

# The problem now is to find the value of beta.
# The additional "swap" condition means that for some time t_s, the coordinates of A become B's and vice versa.
# This geometrical constraint is only satisfied if the common ecliptic longitude lambda_0 is 90 or 270 degrees
# (on the solstitial colure), and if the magnitude of the ecliptic latitude |beta| is equal to epsilon.
# So, |beta| = epsilon.

beta_deg = epsilon_deg
beta_rad = epsilon_rad

# Step 3: Check consistency
# Let's check if this solution (lambda_0=270, beta_A=epsilon, beta_B=-epsilon) works.
# (Note: we could also choose beta_A=-epsilon, beta_B=epsilon, which gives the same distance).
# We must check if star A crosses the equator at t=-3000 years.
# The condition for equator crossing is: sin(lambda_cross) = -tan(beta_A)/tan(epsilon)
# For beta_A = epsilon, tan(beta_A) = tan(epsilon).
# sin(lambda_cross) = -tan(epsilon)/tan(epsilon) = -1.
# This means the crossing longitude must be lambda_cross = 270 degrees.
# Let's check the longitude of star A at t = -3000 years.
# lambda_A(-3000) = lambda_0 - (360/26000) * (-3000)
# lambda_0 is assumed to be 270 degrees.
lambda_0_deg = 270
lambda_A_at_minus_3000 = lambda_0_deg + (360.0 / period) * 3000
# The result should be 270 degrees (mod 360).
# 270 + 360*3/26 = 270 + 41.538... = 311.538...
# This is not 270. Let's try lambda_0 = 90.
lambda_0_deg = 90
# sin(lambda_cross) = -1 still implies lambda_cross = 270.
# lambda_A(-3000) = 90 + 41.538 = 131.538...
# This doesn't seem to work.

# There is a subtle point often missed. Let's re-evaluate the swap.
# The swap condition for stars at (lambda_0, beta) and (lambda_0, -beta) holds true if:
# 1. The time of swap corresponds to a precession of 180 degrees (t_s = 13000 years)
# 2. tan(beta) = tan(epsilon) * sin(lambda_0)
# The equator crossing condition for Star A is:
# sin(lambda_0 + (360/P)*3000) = -tan(beta)/tan(epsilon)
# Substituting the first equation into the second:
# sin(lambda_0 + 41.538...) = -(tan(epsilon)*sin(lambda_0))/tan(epsilon) = -sin(lambda_0)
# sin(lambda_0)cos(41.538) + cos(lambda_0)sin(41.538) = -sin(lambda_0)
# tan(lambda_0)(cos(41.538)+1) = -sin(41.538)
# tan(lambda_0) = -sin(41.538) / (1 + cos(41.538))
# Using the half-angle identity tan(x/2) = sin(x)/(1+cos(x)), this means:
# tan(lambda_0) = -tan(41.538 / 2)
# So, lambda_0 = -20.769 degrees.
precession_angle_A = (360.0/period) * 3000
lambda_0_deg = -precession_angle_A / 2.0
lambda_0_rad = math.radians(lambda_0_deg)

# Now find beta using tan(beta) = tan(epsilon) * sin(lambda_0)
tan_beta = math.tan(epsilon_rad) * math.sin(lambda_0_rad)
beta_rad = math.atan(tan_beta)
beta_deg = math.degrees(beta_rad)

# The angular distance between the stars is 2 * |beta|
angular_distance = 2 * abs(beta_deg)

print("The key timing clue implies that the stars have the same ecliptic longitude and opposite ecliptic latitudes.")
print(f"Let the coordinates be (lambda, beta) and (lambda, -beta).")
print(f"The 'swap' and 'equator crossing' conditions allow us to solve for these.")
print(f"Axial tilt (epsilon) = {epsilon_deg} degrees.")
print(f"Precession angle for Star A's crossing time = (360 / {period}) * {t_A*-1} = {precession_angle_A:.2f} degrees.")
print(f"Solved Ecliptic Longitude (lambda) = {lambda_0_deg:.2f} degrees.")
print(f"Solved Ecliptic Latitude (beta) = {beta_deg:.2f} degrees.")
print(f"The angular distance between the stars is 2 * |beta|.")
print(f"Final Equation: 2 * abs({beta_deg:.2f})")
print(f"Angular Distance = {angular_distance:.1f} degrees.")
