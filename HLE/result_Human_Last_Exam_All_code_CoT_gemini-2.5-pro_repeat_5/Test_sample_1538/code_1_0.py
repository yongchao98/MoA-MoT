import math

def solve_flight_time():
    """
    Calculates the time of flight for a spacecraft in a polar orbit around Mercury.
    """
    # --- Given Parameters ---
    # Convert all units to meters and seconds for consistency
    h_p_km = 218  # Periapsis altitude in km
    h_a_km = 9982 # Apoapsis altitude in km
    R_km = 2440   # Radius of Mercury in km
    mu_m3_s2 = 2.2032e13 # Standard gravitational parameter of Mercury in m^3/s^2
    
    # Convert km to m
    h_p = h_p_km * 1000
    h_a = h_a_km * 1000
    R = R_km * 1000

    print("--- Input Parameters ---")
    print(f"Periapsis altitude (h_p): {h_p_km} km")
    print(f"Apoapsis altitude (h_a): {h_a_km} km")
    print(f"Radius of Mercury (R): {R_km} km")
    print(f"Standard gravitational parameter (mu): {mu_m3_s2:.4e} m^3/s^2\n")

    # --- Step 1: Calculate Orbital Parameters ---
    r_p = R + h_p
    r_a = R + h_a
    a = (r_p + r_a) / 2
    e = (r_a - r_p) / (r_a + r_p)

    print("--- Calculated Orbital Parameters ---")
    print(f"Radius of periapsis (r_p): {int(r_p)} m")
    print(f"Radius of apoapsis (r_a): {int(r_a)} m")
    print(f"Semi-major axis (a): {int(a)} m")
    print(f"Eccentricity (e): {e:.6f}\n")

    # --- Step 2: Determine True Anomaly for Start and End Points ---
    # Argument of periapsis (omega) is the latitude of periapsis for a polar orbit.
    omega_deg = 60.0
    # Start point: North Pole (latitude 90 deg) -> 90 = omega + nu_start -> nu_start = 30 deg
    nu_start_deg = 30.0
    # End point: Equatorial plane (descending node) -> 180 = omega + nu_end -> nu_end = 120 deg
    nu_end_deg = 120.0
    
    nu_start_rad = math.radians(nu_start_deg)
    nu_end_rad = math.radians(nu_end_deg)

    print("--- Flight Path Angles ---")
    print(f"Start true anomaly (nu_start): {nu_start_deg} degrees")
    print(f"End true anomaly (nu_end): {nu_end_deg} degrees\n")

    # --- Step 3: Calculate Time of Flight ---
    
    # Function to calculate mean anomaly from true anomaly
    def get_mean_anomaly(nu_rad, e):
        # Convert true anomaly to eccentric anomaly (E)
        # Note: Using atan2 is more robust but this formula is sufficient for 0 < nu < 180
        E_rad = 2 * math.atan(math.sqrt((1 - e) / (1 + e)) * math.tan(nu_rad / 2))
        # Calculate mean anomaly (M) from eccentric anomaly using Kepler's Equation
        M_rad = E_rad - e * math.sin(E_rad)
        return M_rad

    # Calculate mean anomalies for start and end points
    M_start_rad = get_mean_anomaly(nu_start_rad, e)
    M_end_rad = get_mean_anomaly(nu_end_rad, e)
    
    delta_M_rad = M_end_rad - M_start_rad

    # The final equation for time of flight is: t = sqrt(a^3 / mu) * (M_end - M_start)
    # We will print each component of this equation.
    time_factor = math.sqrt(a**3 / mu_m3_s2)
    time_of_flight_s = time_factor * delta_M_rad

    print("--- Time of Flight Calculation ---")
    print(f"Equation: t = sqrt(a^3 / mu) * (M_end - M_start)")
    print(f"Value of a: {int(a)}")
    print(f"Value of mu: {mu_m3_s2:.4e}")
    print(f"Value of M_start: {M_start_rad:.6f} radians")
    print(f"Value of M_end: {M_end_rad:.6f} radians")
    print(f"Calculated time of flight (unrounded): {time_of_flight_s:.2f} s\n")

    # --- Step 4: Final Answer Formatting ---
    # Round to the nearest 10 seconds
    final_time = int(round(time_of_flight_s / 10.0) * 10)

    print("--- Final Answer ---")
    print(f"The time it took for the spacecraft to travel is {final_time} seconds (rounded to the nearest 10 seconds).")
    
    return final_time

# Execute the function and capture the final answer
final_answer = solve_flight_time()
print(f"<<<{final_answer}>>>")
