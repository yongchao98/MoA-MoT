import math

# Step 1: Define constants and calculate orbital parameters
# Given values
h_p = 218 * 1000        # Periapsis altitude in meters
h_a = 9982 * 1000       # Apoapsis altitude in meters
R_M = 2440 * 1000       # Radius of Mercury in meters
mu = 2.2032e13          # Standard gravitational parameter of Mercury in m^3/s^2

# Calculate radius of periapsis (r_p) and apoapsis (r_a)
r_p = R_M + h_p
r_a = R_M + h_a

# Calculate semi-major axis (a)
a = (r_p + r_a) / 2

# Calculate eccentricity (e)
e = (r_a - r_p) / (r_a + r_p)

# Calculate mean motion (n)
n = math.sqrt(mu / a**3)

print("--- Orbital Parameters ---")
print(f"Semi-major axis (a): {a:.0f} m")
print(f"Eccentricity (e): {e:.6f}")
print(f"Mean motion (n): {n:.6f} rad/s")
print("-" * 26)

# Step 2: Determine true anomalies
# The path is from the North Pole (before periapsis) to the Equatorial plane (after periapsis).
# Based on the geometry (see explanation), this corresponds to:
theta1_deg = -30.0  # True anomaly at start (North Pole)
theta2_deg = 60.0   # True anomaly at end (Equatorial Plane)

theta1_rad = math.radians(theta1_deg)
theta2_rad = math.radians(theta2_deg)

print("\n--- Flight Path Anomalies ---")
print(f"Start True Anomaly (theta_1): {theta1_deg}°")
print(f"End True Anomaly (theta_2): {theta2_deg}°")
print("-" * 26)

# Step 3: Calculate Time of Flight using Kepler's Equation

def get_time_from_periapsis(theta_rad, e, n):
    """Calculates time elapsed since periapsis passage for a given true anomaly."""
    # Convert true anomaly (theta) to eccentric anomaly (E)
    tan_E_over_2 = math.sqrt((1 - e) / (1 + e)) * math.tan(theta_rad / 2)
    E_rad = 2 * math.atan(tan_E_over_2)
    
    # Convert eccentric anomaly (E) to mean anomaly (M)
    M_rad = E_rad - e * math.sin(E_rad)
    
    # Convert mean anomaly (M) to time (t)
    t = M_rad / n
    
    return E_rad, M_rad, t

# Calculate for the start point
E1_rad, M1_rad, t1 = get_time_from_periapsis(theta1_rad, e, n)

# Calculate for the end point
E2_rad, M2_rad, t2 = get_time_from_periapsis(theta2_rad, e, n)

print("\n--- Time Calculation ---")
print("Start Point (theta_1 = -30°):")
print(f"  Eccentric Anomaly (E_1): {math.degrees(E1_rad):.4f}°")
print(f"  Mean Anomaly (M_1): {math.degrees(M1_rad):.4f}°")
print(f"  Time from periapsis (t_1): {t1:.2f} s")

print("\nEnd Point (theta_2 = 60°):")
print(f"  Eccentric Anomaly (E_2): {math.degrees(E2_rad):.4f}°")
print(f"  Mean Anomaly (M_2): {math.degrees(M2_rad):.4f}°")
print(f"  Time from periapsis (t_2): {t2:.2f} s")

# The total time of flight is the difference
time_of_flight = t2 - t1

# Step 4: Final rounding
# Round the final answer to the nearest 10 seconds
final_answer = int(round(time_of_flight / 10.0) * 10)

print("-" * 26)
print(f"\nTotal time of flight = t_2 - t_1 = {t2:.2f} - ({t1:.2f}) = {time_of_flight:.2f} seconds")
print(f"\nFinal Answer (rounded to the nearest 10 seconds): {final_answer}")