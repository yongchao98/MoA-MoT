import math

# Given constants
T = 26000  # years, precession period
t_A_ago = 3000  # years
t_B_future = 10000  # years
epsilon_deg = 23.5  # degrees, axial tilt

# Convert to radians for calculation
epsilon_rad = math.radians(epsilon_deg)

# In the derivation, we found a relationship between the angle phi and the timing data.
# The L-motion in 10000 years must be less than the shortest angular distance between equator crossings (180 - 2*phi).
# The boundary case implies they are equal.
# The angular motion of L in 10000 years is:
delta_L_B_deg = 360.0 * t_B_future / T

# From the boundary condition 2*phi = 180 - delta_L_B_deg
phi_deg = (180.0 - delta_L_B_deg) / 2.0
phi_rad = math.radians(phi_deg)

# The ecliptic latitude beta is related to phi and epsilon by tan(beta) = sin(phi) * tan(epsilon)
tan_beta = math.sin(phi_rad) * math.tan(epsilon_rad)
beta_rad = math.atan(tan_beta)
beta_deg = math.degrees(beta_rad)

# The angular distance 'd' between the stars is 2*beta
angular_distance_deg = 2 * beta_deg

print(f"Period of precession T = {T} years")
print(f"Axial tilt epsilon = {epsilon_deg} degrees")
print(f"Time for Star A since last equator crossing = {t_A_ago} years")
print(f"Time for Star B until first equator crossing = {t_B_future} years")
print(f"Calculation Steps:")
print(f"1. Angular motion of the celestial grid in 10000 years: delta_L_B = 360 * {t_B_future} / {T} = {delta_L_B_deg:.2f} degrees")
print(f"2. From the 'first crossing' boundary condition, we find 2*phi = 180 - {delta_L_B_deg:.2f} = {2*phi_deg:.2f} degrees")
print(f"3. This gives the angle phi = {phi_deg:.2f} degrees")
print(f"4. The ecliptic latitude beta is found using tan(beta) = tan(epsilon) * sin(phi)")
print(f"   beta = arctan(tan({epsilon_deg}) * sin({phi_deg:.2f})) = {beta_deg:.2f} degrees")
print(f"5. The angular distance between the stars is d = 2 * beta")
print(f"   Final Equation: d = 2 * {beta_deg:.2f}")
print(f"Resulting angular distance: {angular_distance_deg:.2f} degrees")
