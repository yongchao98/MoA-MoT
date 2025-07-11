import math

def solve_orbital_time():
    """
    Calculates the time for a spacecraft to travel from the north pole
    through periapsis to the equatorial plane of Mercury.
    """
    # Given constants
    h_p_km = 218
    h_a_km = 9982
    R_km = 2440
    mu = 2.2032e13  # m^3/s^2

    # Convert all units to SI (meters)
    h_p = h_p_km * 1000
    h_a = h_a_km * 1000
    R = R_km * 1000

    # Step 1: Calculate orbital parameters
    r_p = R + h_p
    r_a = R + h_a
    a = (r_p + r_a) / 2
    e = (r_a - r_p) / (r_a + r_p)

    # Step 2: Define start and end true anomalies based on the problem geometry
    # Start point (North Pole, 90N) corresponds to nu_1 = -30 degrees
    # End point (Equator, 0N) corresponds to nu_2 = 60 degrees
    nu1_deg = -30.0
    nu2_deg = 60.0

    # Step 3: Function to calculate time from periapsis using Kepler's Equation
    def get_time_from_periapsis(nu_deg, a_val, e_val, mu_val):
        """Calculates time from periapsis for a given true anomaly."""
        nu_rad = math.radians(nu_deg)
        
        # Calculate eccentric anomaly (E) from true anomaly (nu)
        E_half = math.atan(math.sqrt((1 - e_val) / (1 + e_val)) * math.tan(nu_rad / 2))
        E = 2 * E_half
        
        # Calculate mean anomaly (M) from eccentric anomaly (E)
        M = E - e_val * math.sin(E)
        
        # Calculate time from periapsis (t) from mean anomaly (M)
        time = M * math.sqrt(a_val**3 / mu_val)
        return time

    # Calculate time relative to periapsis for the start and end points
    t1 = get_time_from_periapsis(nu1_deg, a, e, mu)
    t2 = get_time_from_periapsis(nu2_deg, a, e, mu)

    # Step 4: Calculate total time of flight and round
    total_time = t2 - t1
    rounded_time = int(round(total_time / 10) * 10)

    # Print the detailed calculation including the final equation
    print("--- Calculation Details ---")
    print(f"Semi-major axis (a): {a:.0f} m")
    print(f"Eccentricity (e): {e:.6f}")
    print("-" * 25)
    print(f"Time to reach start position (nu={nu1_deg}°): {t1:.2f} s")
    print(f"Time to reach end position (nu={nu2_deg}°): {t2:.2f} s")
    print("-" * 25)
    print("Final Time of Flight Equation:")
    print(f"Total Time = (Time to reach {nu2_deg}°) - (Time to reach {nu1_deg}°)")
    print(f"Total Time = {t2:.2f} s - ({t1:.2f} s)")
    print(f"Total Time = {total_time:.2f} s")
    print("-" * 25)
    print(f"Final Answer (rounded to the nearest 10 seconds): {rounded_time}")

solve_orbital_time()
<<<1270>>>