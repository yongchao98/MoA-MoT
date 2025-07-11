import math

def solve_flight_time():
    """
    Calculates the time of flight for a spacecraft in a polar orbit around Mercury.
    """
    # --- Step 1: Define constants and convert to SI units ---
    # Given values
    h_p_km = 218
    h_a_km = 9982
    R_M_km = 2440
    mu_M_m3_s2 = 2.2032e13

    # Convert km to meters
    m_per_km = 1000
    h_p = h_p_km * m_per_km
    h_a = h_a_km * m_per_km
    R_M = R_M_km * m_per_km
    mu_M = mu_M_m3_s2

    # --- Step 2: Calculate orbital parameters ---
    # Radii of periapsis and apoapsis
    r_p = R_M + h_p
    r_a = R_M + h_a

    # Semi-major axis (a)
    a = (r_p + r_a) / 2.0

    # Eccentricity (e)
    e = (r_a - r_p) / (r_a + r_p)

    # Mean motion (n) in rad/s
    n = math.sqrt(mu_M / a**3)

    # --- Step 3: Define true anomalies for the flight path ---
    # The path is from the North Pole (before periapsis) to the equatorial plane (after periapsis).
    # For this polar orbit, this corresponds to true anomalies from -30 to +60 degrees.
    nu1_rad = math.radians(-30)
    nu2_rad = math.radians(60)

    # --- Step 4: Convert true anomalies to mean anomalies ---
    def true_to_mean(nu_rad, e):
        """Converts true anomaly to mean anomaly."""
        # First, find the eccentric anomaly (E)
        cos_E = (e + math.cos(nu_rad)) / (1 + e * math.cos(nu_rad))
        sin_E = (math.sqrt(1 - e**2) * math.sin(nu_rad)) / (1 + e * math.cos(nu_rad))
        E_rad = math.atan2(sin_E, cos_E)
        
        # Then, find the mean anomaly (M) using Kepler's Equation
        M_rad = E_rad - e * math.sin(E_rad)
        return M_rad

    M1_rad = true_to_mean(nu1_rad, e)
    M2_rad = true_to_mean(nu2_rad, e)

    # --- Step 5: Calculate the time of flight ---
    # The change in mean anomaly divided by mean motion
    delta_M = M2_rad - M1_rad
    delta_t = delta_M / n

    # --- Step 6: Print the final calculation and round the answer ---
    print(f"The final calculation for time of flight (t) is (M2 - M1) / n:")
    print(f"t = ({M2_rad:.6f} rad - ({M1_rad:.6f} rad)) / {n:.9f} rad/s")
    print(f"t = {delta_t:.2f} s")

    # Round the result to the nearest 10 seconds
    rounded_time = int(round(delta_t / 10) * 10)
    print(f"\nFinal Answer (rounded to the nearest 10 seconds): {rounded_time}")
    
    return rounded_time

# Execute the function and capture the final answer for the 'answer block'
final_answer = solve_flight_time()
print(f"<<<{final_answer}>>>")
