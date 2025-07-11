import math

# Step 1: Define constants and initial values in SI units
R_M = 2440 * 1000  # Radius of Mercury in meters
h_p = 218 * 1000   # Periapsis altitude in meters
h_a = 9982 * 1000  # Apoapsis altitude in meters
mu = 2.2032e13     # Standard gravitational parameter of Mercury in m^3/s^2

# Step 2: Calculate orbital parameters
r_p = R_M + h_p  # Periapsis radius
r_a = R_M + h_a  # Apoapsis radius
a = (r_p + r_a) / 2  # Semi-major axis
e = (r_a - r_p) / (r_a + r_p)  # Eccentricity

# Step 3: Define start and end true anomalies based on the path
# Path: North Pole -> Periapsis -> Equator
# This corresponds to nu_start = -30 deg and nu_end = 60 deg
nu_start_deg = -30.0
nu_end_deg = 60.0

# Convert true anomalies to radians for calculations
nu_start_rad = math.radians(nu_start_deg)
nu_end_rad = math.radians(nu_end_deg)

# Step 4: Calculate time of flight using Kepler's Equation

def get_mean_anomaly(true_anomaly_rad, eccentricity):
    """Calculates mean anomaly from true anomaly and eccentricity."""
    # Calculate eccentric anomaly (E)
    E_half_tan_val = math.sqrt((1 - eccentricity) / (1 + eccentricity)) * math.tan(true_anomaly_rad / 2)
    E_rad = 2 * math.atan(E_half_tan_val)
    
    # Calculate mean anomaly (M)
    M_rad = E_rad - eccentricity * math.sin(E_rad)
    return M_rad

# Calculate mean anomalies for start and end points
M_start_rad = get_mean_anomaly(nu_start_rad, e)
M_end_rad = get_mean_anomaly(nu_end_rad, e)

# Calculate mean motion (n)
n = math.sqrt(mu / (a**3))

# Calculate time of flight (delta_t)
delta_M = M_end_rad - M_start_rad
delta_t = delta_M / n

# Step 5: Round the result to the nearest 10 seconds
rounded_time = int(round(delta_t / 10) * 10)

# --- Final Output ---
print("This script calculates the time of flight for a spacecraft in orbit around Mercury.")
print("\n--- Orbital Parameters ---")
print(f"Periapsis radius (r_p): {r_p / 1000:.1f} km")
print(f"Apoapsis radius (r_a): {r_a / 1000:.1f} km")
print(f"Semi-major axis (a): {a / 1000:.1f} km")
print(f"Eccentricity (e): {e:.6f}")

print("\n--- Time of Flight Calculation ---")
print(f"Start true anomaly (ν_start): {nu_start_deg}° -> Start Mean Anomaly (M_start): {M_start_rad:.6f} rad")
print(f"End true anomaly (ν_end): {nu_end_deg}° -> End Mean Anomaly (M_end): {M_end_rad:.6f} rad")
print(f"Mean Anomaly Difference (ΔM) = M_end - M_start = {M_end_rad:.6f} - ({M_start_rad:.6f}) = {delta_M:.6f} rad")
print(f"Mean Motion (n) = sqrt(μ / a^3) = sqrt({mu:.4e} / {a**3:.4e}) = {n:.8f} rad/s")
print("\nFinal Equation for Time:")
print(f"Time of Flight (Δt) = ΔM / n = {delta_M:.6f} rad / {n:.8f} rad/s = {delta_t:.1f} seconds")

print(f"\nThe calculated time is {delta_t:.1f} seconds.")
print(f"The time rounded to the nearest 10 seconds is: {rounded_time}")
<<<12690>>>