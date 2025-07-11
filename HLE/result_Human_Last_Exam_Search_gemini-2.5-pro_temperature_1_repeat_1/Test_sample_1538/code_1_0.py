import math

# Given constants
h_p = 218 * 1000  # periapsis altitude in meters
h_a = 9982 * 1000 # apoapsis altitude in meters
R = 2440 * 1000   # radius of Mercury in meters
mu = 2.2032e13    # standard gravitational parameter of Mercury in m^3/s^2

# 1. Calculate orbital parameters
r_p = R + h_p  # periapsis radius
r_a = R + h_a  # apoapsis radius

a = (r_p + r_a) / 2  # semi-major axis
e = (r_a - r_p) / (r_a + r_p) # eccentricity
n = math.sqrt(mu / a**3) # mean motion

# 2. Determine true anomalies
# From the problem statement and geometric analysis:
# The path is from North Pole (90 deg lat) -> Periapsis (60 deg lat) -> Equator (0 deg lat)
# This requires an argument of periapsis (omega) of 120 degrees.
# This gives:
# True anomaly at North Pole (nu1) = -30 degrees
# True anomaly at Equator (nu2) = 60 degrees

nu1_deg = -30.0
nu2_deg = 60.0

nu1_rad = math.radians(nu1_deg)
nu2_rad = math.radians(nu2_deg)

# 3. Calculate time of flight
def calculate_time_from_periapsis(nu_rad, e, n):
    """Calculates the time of flight from periapsis to a given true anomaly."""
    # Convert true anomaly to eccentric anomaly (E)
    # tan(E/2) = sqrt((1-e)/(1+e)) * tan(nu/2)
    tan_E_over_2 = math.sqrt((1 - e) / (1 + e)) * math.tan(nu_rad / 2)
    E_rad = 2 * math.atan(tan_E_over_2)

    # Calculate mean anomaly (M) using Kepler's Equation
    # M = E - e * sin(E)
    M_rad = E_rad - e * math.sin(E_rad)

    # Calculate time (t)
    # M = n * t => t = M / n
    t = M_rad / n
    return t

# Time at the North Pole (initial point)
t1 = calculate_time_from_periapsis(nu1_rad, e, n)

# Time at the Equator (final point)
t2 = calculate_time_from_periapsis(nu2_rad, e, n)

# Total time of flight is the difference
total_time = t2 - t1

# Round to the nearest 10 seconds
rounded_time = int(round(total_time / 10.0) * 10)

# print(f"Periapsis radius (r_p): {r_p} m")
# print(f"Apoapsis radius (r_a): {r_a} m")
# print(f"Semi-major axis (a): {a} m")
# print(f"Eccentricity (e): {e}")
# print(f"Mean motion (n): {n} rad/s")
# print(f"True anomaly at start (nu1): {nu1_deg} deg")
# print(f"True anomaly at end (nu2): {nu2_deg} deg")
# print(f"Time at start (t1): {t1} s")
# print(f"Time at end (t2): {t2} s")
# print(f"Total time of flight: {total_time} s")
# print(f"Rounded time: {rounded_time} s")

# Final Answer
# print(f"<<<{rounded_time}>>>")