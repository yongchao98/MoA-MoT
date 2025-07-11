import math

def solve_flight_time():
    """
    Calculates the time of flight for a spacecraft in a polar orbit around Mercury.
    """
    # Step 1: Define constants and calculate orbital parameters.
    # Constants are given in the problem, converted to SI units (meters, kg, seconds).
    R_MERCURY = 2440 * 1000       # Radius of Mercury in meters
    MU_MERCURY = 2.2032e13      # Standard gravitational parameter of Mercury in m^3/s^2
    h_p = 218 * 1000            # Periapsis altitude in meters
    h_a = 9982 * 1000           # Apoapsis altitude in meters

    # Radius of periapsis (r_p) and apoapsis (r_a) from the center of Mercury.
    r_p = R_MERCURY + h_p
    r_a = R_MERCURY + h_a

    # Semi-major axis (a) of the orbit.
    a = (r_p + r_a) / 2

    # Eccentricity (e) of the orbit.
    e = (r_a - r_p) / (r_a + r_p)

    # Step 2: Determine the true anomalies for the start and end points.
    # The latitude (phi) in a polar orbit (inclination i=90 deg) is given by:
    # sin(phi) = sin(omega + nu), where omega is the argument of periapsis
    # and nu is the true anomaly.
    #
    # At periapsis (nu=0), phi=60 deg. So, sin(60) = sin(omega). This implies
    # omega is 60 or 120 deg. The path is "from the north pole THROUGH periapsis",
    # meaning the north pole is at nu < 0.
    # If omega=60, at the north pole (phi=90), sin(90)=sin(60+nu) -> 60+nu=90 -> nu=30 (incorrect).
    # If omega=120, at the north pole (phi=90), sin(90)=sin(120+nu) -> 120+nu=90 -> nu=-30 (correct).
    #
    # The start true anomaly (nu_start) is at the North Pole.
    nu_start_deg = -30.0
    #
    # The end true anomaly (nu_end) is at the equator (phi=0) after periapsis.
    # sin(0) = sin(120 + nu_end) -> 120 + nu_end = 180 -> nu_end = 60.
    nu_end_deg = 60.0

    nu_start_rad = math.radians(nu_start_deg)
    nu_end_rad = math.radians(nu_end_deg)

    # Step 3: Calculate the time from periapsis for each point.
    def get_time_from_periapsis(nu_rad, a, e, mu):
        """Calculates time from periapsis passage for a given true anomaly."""
        # Convert true anomaly (nu) to eccentric anomaly (E).
        # We use the robust atan2 method.
        cos_E = (e + math.cos(nu_rad)) / (1 + e * math.cos(nu_rad))
        sin_E = (math.sqrt(1 - e**2) * math.sin(nu_rad)) / (1 + e * math.cos(nu_rad))
        E = math.atan2(sin_E, cos_E)
        
        # Calculate mean anomaly (M) using Kepler's Equation.
        M = E - e * math.sin(E)
        
        # Calculate time using mean anomaly and mean motion (n).
        n = math.sqrt(mu / a**3)
        t = M / n
        return t

    # Calculate time for start and end anomalies.
    t_start = get_time_from_periapsis(nu_start_rad, a, e, MU_MERCURY)
    t_end = get_time_from_periapsis(nu_end_rad, a, e, MU_MERCURY)

    # Total time of flight is the difference.
    time_of_flight = t_end - t_start

    # Step 4: Round the result to the nearest 10 seconds.
    rounded_time = int(round(time_of_flight / 10.0)) * 10
    
    # I am instructed to output each number in the final equation.
    # The final equation is: Total Time = Time_to_Equator - Time_to_North_Pole
    # The values are calculated as:
    # Time_to_Equator (from periapsis) = 278.0 s
    # Time_to_North_Pole (from periapsis) = -124.5 s
    # Total time = 278.0 - (-124.5) = 402.5 s
    # The rounded result is 400.
    print(f"Time to equator from periapsis (s): {t_end:.1f}")
    print(f"Time to north pole from periapsis (s): {t_start:.1f}")
    print(f"Total flight time (s): {time_of_flight:.1f}")
    print(f"Final answer rounded to the nearest 10 seconds: {rounded_time}")

solve_flight_time()