import math

def solve_orbital_time():
    """
    Calculates the time for a spacecraft to travel from the North Pole
    through periapsis to the equatorial plane around Mercury.
    """
    # Step 1: Define constants and orbital parameters
    R_M = 2440e3  # Radius of Mercury in meters
    mu = 2.2032e13  # Standard gravitational parameter of Mercury in m^3/s^2
    h_p = 218e3  # Periapsis altitude in meters
    h_a = 9982e3 # Apoapsis altitude in meters

    # Radii from the center of Mercury
    r_p = R_M + h_p
    r_a = R_M + h_a

    # Semi-major axis (a) and eccentricity (e)
    a = (r_p + r_a) / 2
    e = (r_a - r_p) / (r_a + r_p)

    # Step 2: This step is the logical derivation of the true anomalies,
    # which are found to be 30 and -60 degrees. By symmetry, we calculate the time
    # for positive angles 30 and 60 degrees.
    nu1_deg = 30
    nu2_deg = 60

    def get_time_from_periapsis(true_anomaly_deg):
        """Calculates time from periapsis passage to a given true anomaly."""
        nu_rad = math.radians(true_anomaly_deg)

        # Calculate eccentric anomaly (E) from true anomaly (nu)
        # Using the robust atan2 form to avoid quadrant ambiguity
        cos_E = (e + math.cos(nu_rad)) / (1 + e * math.cos(nu_rad))
        sin_E = (math.sqrt(1 - e**2) * math.sin(nu_rad)) / (1 + e * math.cos(nu_rad))
        E_rad = math.atan2(sin_E, cos_E)
        
        # Calculate mean anomaly (M) from eccentric anomaly (E)
        # M = E - e*sin(E)
        M = E_rad - e * math.sin(E_rad)

        # Calculate time (t) from mean anomaly
        # M = n*t, where n = sqrt(mu/a^3) is the mean motion
        t = M * math.sqrt(a**3 / mu)
        return t

    # Step 3: Calculate the two time intervals
    # Time from North Pole (nu=30 deg) to periapsis
    time1 = get_time_from_periapsis(nu1_deg)
    # Time from periapsis to Equator (nu=-60 deg)
    time2 = get_time_from_periapsis(nu2_deg)
    
    # Total time is the sum of the two intervals
    total_time = time1 + time2

    # Step 4: Round to the nearest 10 seconds
    rounded_time = int(round(total_time / 10.0) * 10)

    # Print the breakdown of the calculation as requested
    print("This problem asks for the time it took for the spacecraft to travel from the north pole through periapsis to the equatorial plane.")
    print(f"The first part of the journey is from the north pole (true anomaly = {nu1_deg} degrees) to periapsis.")
    print(f"The second part is from periapsis to the equatorial plane (true anomaly = -{nu2_deg} degrees).")
    print(f"\nUsing Kepler's equation, we calculate the time for each segment:")
    print(f"Time from periapsis to true anomaly {nu1_deg} degrees = {int(round(time1))} seconds")
    print(f"Time from periapsis to true anomaly {nu2_deg} degrees = {int(round(time2))} seconds")
    print(f"\nThe final equation for the total time is:")
    print(f"{int(round(time1))} + {int(round(time2))} = {int(round(total_time))} seconds")
    print(f"\nRounding to the nearest 10 seconds, the final answer is {rounded_time}.")

    # Final answer in the specified format
    print(f"\n<<<{rounded_time}>>>")

solve_orbital_time()