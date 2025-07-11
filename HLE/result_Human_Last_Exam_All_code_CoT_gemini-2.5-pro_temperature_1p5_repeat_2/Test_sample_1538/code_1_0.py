import math

def solve_orbital_time():
    """
    Calculates the time for a spacecraft to travel from the North Pole, through
    periapsis, to the equatorial plane of Mercury.
    """
    # Step 1: Define constants in SI units
    R_mercury = 2440 * 1000  # Radius of Mercury in meters
    mu_mercury = 2.2032e13   # Standard gravitational parameter of Mercury in m^3/s^2
    h_p = 218 * 1000         # Periapsis altitude in meters
    h_a = 9982 * 1000        # Apoapsis altitude in meters

    # Step 2: Calculate orbital parameters
    r_p = R_mercury + h_p
    r_a = R_mercury + h_a
    
    a = (r_p + r_a) / 2
    e = (r_a - r_p) / (r_a + r_p)
    n = math.sqrt(mu_mercury / a**3)

    # Step 3: Define true anomalies for the start and end points
    # Periapsis is at 60 deg N latitude. For a polar orbit, lat = omega + nu.
    # At periapsis, nu=0, so the argument of periapsis omega = 60 deg.
    # Start point: North Pole (90 deg N) => 90 = 60 + nu_1 => nu_1 = 30 deg.
    # End point: Equatorial plane (0 deg) after periapsis => 0 = 60 + nu_2 => nu_2 = -60 deg.
    nu1_deg = 30.0
    nu2_deg = -60.0

    def get_time_from_periapsis(nu_deg, e, n):
        """
        Calculates the time elapsed since periapsis for a given true anomaly using Kepler's Equation.
        The time value is always positive.
        """
        nu_rad = math.radians(nu_deg)
        
        # Calculate eccentric anomaly (E)
        term = math.sqrt((1 - e) / (1 + e)) * math.tan(nu_rad / 2)
        E_rad = 2 * math.atan(term)
        
        # Calculate mean anomaly (M)
        M_rad = E_rad - e * math.sin(E_rad)
        
        # Calculate time (t = M/n)
        # We take the absolute value as we are calculating time segments.
        time = abs(M_rad / n)
        return time

    # Step 4: Calculate the time for each segment of the journey
    # Time from North Pole (nu=30) to periapsis (nu=0)
    time_segment1 = get_time_from_periapsis(nu1_deg, e, n)
    
    # Time from periapsis (nu=0) to equatorial plane (nu=-60)
    time_segment2 = get_time_from_periapsis(nu2_deg, e, n)

    # Total time is the sum of the two segments
    total_time = time_segment1 + time_segment2
    
    # Round the final answer to the nearest 10 seconds
    rounded_time = int(round(total_time / 10.0)) * 10

    # Step 5: Print the results
    print("The final equation for the total time is:")
    print(f"Total Time = Time(ν=0 -> ν={nu1_deg}°) + Time(ν=0 -> ν={abs(nu2_deg)}°)")
    print(f"Total Time = {time_segment1} + {time_segment2}")
    print(f"Total Time = {total_time} seconds")
    print("\nFinal Answer:")
    print(f"The total time rounded to the nearest 10 seconds is: {rounded_time}")

solve_orbital_time()
<<<1270>>>