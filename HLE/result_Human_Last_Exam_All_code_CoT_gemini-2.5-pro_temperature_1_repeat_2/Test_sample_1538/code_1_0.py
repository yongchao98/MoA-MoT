import math

def solve_flight_time():
    """
    Calculates the time of flight for a spacecraft in a polar orbit around Mercury.
    """
    # 1. Given constants and conversion to SI units (meters, seconds)
    h_p_km = 218
    h_a_km = 9982
    R_mercury_km = 2440
    mu_m3_s2 = 2.2032e13

    h_p = h_p_km * 1000  # Periapsis altitude in meters
    h_a = h_a_km * 1000  # Apoapsis altitude in meters
    R_mercury = R_mercury_km * 1000 # Mercury radius in meters

    # 2. Calculate orbital parameters
    r_p = R_mercury + h_p  # Periapsis radius
    r_a = R_mercury + h_a  # Apoapsis radius
    
    a = (r_p + r_a) / 2  # Semi-major axis
    e = (r_a - r_p) / (r_a + r_p)  # Eccentricity
    
    # Mean motion (n) in rad/s
    n = math.sqrt(mu_m3_s2 / a**3)

    # 3. Define true anomalies for the start and end of the trajectory
    # From the problem analysis:
    # Start at North Pole (latitude 90 deg) -> nu_1 = -30 deg
    # End at Equatorial Plane (latitude 0 deg) -> nu_2 = 60 deg
    nu1_deg = -30.0
    nu2_deg = 60.0

    # 4. Calculate time of flight using Kepler's Equation
    def get_time_from_periapsis(nu_deg, e, n):
        """Calculates time from periapsis passage for a given true anomaly."""
        nu_rad = math.radians(nu_deg)
        
        # Calculate eccentric anomaly (E)
        tan_E_half = math.sqrt((1 - e) / (1 + e)) * math.tan(nu_rad / 2)
        E_rad = 2 * math.atan(tan_E_half)
        
        # Calculate mean anomaly (M) using Kepler's Equation
        M_rad = E_rad - e * math.sin(E_rad)
        
        # Calculate time from periapsis (t = M/n)
        time = M_rad / n
        return time

    # Calculate time for start and end points
    t1 = get_time_from_periapsis(nu1_deg, e, n)
    t2 = get_time_from_periapsis(nu2_deg, e, n)

    # Total time of flight is the difference
    time_of_flight = t2 - t1

    # 5. Round the result to the nearest 10 seconds
    rounded_time_of_flight = int(round(time_of_flight / 10) * 10)

    # Output the steps and the final answer
    print(f"Orbital Parameters:")
    print(f"  Periapsis radius (r_p) = {r_p:.0f} m")
    print(f"  Apoapsis radius (r_a) = {r_a:.0f} m")
    print(f"  Semi-major axis (a) = {a:.0f} m")
    print(f"  Eccentricity (e) = {e:.6f}")
    print(f"  Mean motion (n) = {n:.6e} rad/s\n")
    
    print(f"Trajectory Calculation:")
    print(f"  Start True Anomaly (nu_1) = {nu1_deg} degrees")
    print(f"  End True Anomaly (nu_2) = {nu2_deg} degrees\n")

    print(f"Time from Periapsis Calculation:")
    print(f"  Time to reach start point from periapsis (t1) = {t1:.1f} s")
    print(f"  Time to reach end point from periapsis (t2) = {t2:.1f} s\n")

    print(f"Total Time of Flight Equation:")
    print(f"  Time = t2 - t1")
    print(f"  Time = {t2:.1f} - ({t1:.1f})")
    print(f"  Time = {time_of_flight:.1f} s\n")

    print(f"Final Answer:")
    print(f"Time rounded to the nearest 10 seconds = {rounded_time_of_flight} s")
    
    # Final answer in the required format
    print(f"\n<<<{rounded_time_of_flight}>>>")

# Execute the function
solve_flight_time()