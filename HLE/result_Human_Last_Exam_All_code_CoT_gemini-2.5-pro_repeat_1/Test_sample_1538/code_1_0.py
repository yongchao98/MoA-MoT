import numpy as np

def solve_spacecraft_time():
    """
    Calculates the time it took for a spacecraft to travel from the north pole,
    through periapsis, to the equatorial plane of Mercury.
    """
    # Step 1: Define constants and convert units
    h_p = 218 * 1000  # Periapsis altitude in meters
    h_a = 9982 * 1000  # Apoapsis altitude in meters
    R = 2440 * 1000  # Radius of Mercury in meters
    mu = 2.2032e13  # Standard gravitational parameter of Mercury in m^3/s^2

    # Step 2: Calculate orbital parameters
    r_p = R + h_p
    r_a = R + h_a
    a = (r_p + r_a) / 2
    e = (r_a - r_p) / (r_a + r_p)

    print("--- Orbital Parameters ---")
    print(f"Periapsis radius (r_p): {r_p / 1000:.1f} km")
    print(f"Apoapsis radius (r_a): {r_a / 1000:.1f} km")
    print(f"Semi-major axis (a): {a / 1000:.1f} km")
    print(f"Eccentricity (e): {e:.5f}")
    print("-" * 26)

    # Step 3: Define true anomalies for the start and end points
    # As derived in the plan, the path starts at theta = -30 deg and ends at theta = 60 deg.
    theta_start_deg = -30
    theta_end_deg = 60
    theta_start_rad = np.deg2rad(theta_start_deg)
    theta_end_rad = np.deg2rad(theta_end_deg)
    
    print("--- Flight Path Angles ---")
    print(f"Start True Anomaly (North Pole): {theta_start_deg} degrees")
    print(f"End True Anomaly (Equator): {theta_end_deg} degrees")
    print("-" * 26)

    # Step 4: Calculate time from periapsis for each point
    def get_time_from_periapsis(theta_rad, a, e, mu):
        """Calculates time from periapsis for a given true anomaly."""
        # Calculate eccentric anomaly (E)
        tan_E_half = np.sqrt((1 - e) / (1 + e)) * np.tan(theta_rad / 2)
        E = 2 * np.arctan(tan_E_half)
        
        # Calculate mean anomaly (M)
        M = E - e * np.sin(E)
        
        # Calculate mean motion (n)
        n = np.sqrt(mu / a**3)
        
        # Calculate time from periapsis (t)
        t = M / n
        return t, E, M

    t_start, E_start, M_start = get_time_from_periapsis(theta_start_rad, a, e, mu)
    t_end, E_end, M_end = get_time_from_periapsis(theta_end_rad, a, e, mu)
    
    print("--- Calculation for Start Point ---")
    print(f"Eccentric Anomaly (E_start): {np.rad2deg(E_start):.2f} degrees")
    print(f"Mean Anomaly (M_start): {np.rad2deg(M_start):.2f} degrees")
    print(f"Time from Periapsis (t_start): {t_start:.2f} s")
    print("-" * 33)

    print("--- Calculation for End Point ---")
    print(f"Eccentric Anomaly (E_end): {np.rad2deg(E_end):.2f} degrees")
    print(f"Mean Anomaly (M_end): {np.rad2deg(M_end):.2f} degrees")
    print(f"Time from Periapsis (t_end): {t_end:.2f} s")
    print("-" * 31)

    # Calculate total time of flight
    total_time = t_end - t_start
    
    print("--- Final Time Calculation ---")
    print(f"Total time = t_end - t_start")
    print(f"Total time = {t_end:.2f} - ({t_start:.2f})")
    print(f"Total time = {total_time:.2f} s")
    
    # Step 5: Round the result to the nearest 10 seconds
    rounded_time = int(round(total_time / 10) * 10)
    
    print("\nFinal Answer (rounded to the nearest 10 seconds):")
    print(rounded_time)
    
    return rounded_time

# Execute the function to get the final answer for the submission format
final_answer = solve_spacecraft_time()
print(f"<<<{final_answer}>>>")
