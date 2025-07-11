import math

def solve_flight_time():
    """
    Calculates the time for a spacecraft to travel from the North Pole to the
    equator in a specified orbit around Mercury.
    """
    # 1. Define constants in SI units (meters, seconds)
    R_MERCURY = 2440 * 1000  # Radius of Mercury in meters
    MU_MERCURY = 2.2032e13  # Standard gravitational parameter of Mercury in m^3/s^2
    H_P = 218 * 1000  # Periapsis altitude in meters
    H_A = 9982 * 1000  # Apoapsis altitude in meters

    # 2. Calculate key orbital parameters
    # Radii of periapsis and apoapsis
    r_p = R_MERCURY + H_P
    r_a = R_MERCURY + H_A

    # Semi-major axis (a)
    a = (r_p + r_a) / 2

    # Eccentricity (e)
    e = (r_a - r_p) / (r_a + r_p)

    # Mean motion (n) in rad/s
    n = math.sqrt(MU_MERCURY / a**3)
    
    # 3. Define the true anomalies for the start and end points of the travel segments.
    # The path is from the North Pole (latitude 90 deg) to the Equator (latitude 0 deg)
    # through periapsis (at latitude 60 deg).
    # For a polar orbit, latitude phi = omega + nu, where omega is argument of periapsis.
    # At periapsis (nu=0), phi=60 deg, so omega=60 deg.
    # Start (North Pole, phi=90): 90 = 60 + nu_start => nu_start = 30 deg.
    # End (Equator, phi=0): 0 = 60 + nu_end => nu_end = -60 deg.
    # The total travel is composed of two segments:
    # Segment 1: North Pole (nu=30 deg) to Periapsis (nu=0 deg)
    # Segment 2: Periapsis (nu=0 deg) to Equator (nu=-60 deg)
    
    nu1_deg = 30.0  # True anomaly at North Pole
    nu2_deg = 60.0  # Absolute true anomaly at Equator crossing

    def get_time_from_periapsis(nu_deg, e, n):
        """Calculates time from periapsis for a given true anomaly."""
        nu_rad = math.radians(nu_deg)
        
        # Calculate eccentric anomaly (E)
        tan_E_over_2 = math.sqrt((1 - e) / (1 + e)) * math.tan(nu_rad / 2)
        E_rad = 2 * math.atan(tan_E_over_2)
        
        # Calculate mean anomaly (M)
        M_rad = E_rad - e * math.sin(E_rad)
        
        # Calculate time from periapsis passage
        t = M_rad / n
        return t

    # 4. Calculate time for each segment
    # Time from North Pole to Periapsis (nu=30 deg to nu=0 deg), which has the same duration
    # as traveling from periapsis to nu=30 deg.
    time1 = get_time_from_periapsis(nu1_deg, e, n)
    
    # Time from Periapsis to Equator (nu=0 deg to nu=-60 deg), which has the same duration
    # as traveling from periapsis to nu=60 deg.
    time2 = get_time_from_periapsis(nu2_deg, e, n)

    # Total time is the sum of the two segments
    total_time_s = time1 + time2
    
    # 5. Round the final answer to the nearest 10 seconds
    rounded_time = int(round(total_time_s / 10.0) * 10)

    # Print the results
    print("Orbital Parameters:")
    print(f"  Periapsis radius (r_p): {r_p / 1000:.1f} km")
    print(f"  Apoapsis radius (r_a): {r_a / 1000:.1f} km")
    print(f"  Semi-major axis (a): {a / 1000:.1f} km")
    print(f"  Eccentricity (e): {e:.6f}")
    print(f"  Mean motion (n): {n:.8f} rad/s\n")
    
    print("Time of Flight Calculation:")
    print(f"  True anomaly at North Pole: {nu1_deg}°")
    print(f"  True anomaly at Equator: -{nu2_deg}°")
    print(f"  Time from North Pole to Periapsis = {time1:.2f} s")
    print(f"  Time from Periapsis to Equator = {time2:.2f} s")
    print(f"  Total time = {time1:.2f} + {time2:.2f} = {total_time_s:.2f} s\n")
    
    print(f"Final Answer (rounded to nearest 10 seconds): {rounded_time}")

solve_flight_time()
<<<1270>>>