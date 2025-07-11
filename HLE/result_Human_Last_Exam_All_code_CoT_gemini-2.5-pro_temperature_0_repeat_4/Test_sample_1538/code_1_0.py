import math

def solve_orbital_time():
    """
    Calculates the time for a spacecraft to travel from the North Pole,
    through periapsis, to the equatorial plane of Mercury.
    """
    # Step 1: Define constants and convert units
    R_mercury_km = 2440
    mu_m3_s2 = 2.2032e13
    h_p_km = 218
    h_a_km = 9982

    # Convert all lengths to meters
    R = R_mercury_km * 1000
    h_p = h_p_km * 1000
    h_a = h_a_km * 1000

    # Step 2: Calculate orbital parameters
    r_p = R + h_p
    r_a = R + h_a
    a = (r_p + r_a) / 2
    e = (r_a - r_p) / (r_a + r_p)

    # Step 3: Determine true anomalies for the start and end points
    # Based on the geometry (polar orbit, periapsis at 60N, path through 90N -> periapsis -> 0N),
    # the argument of periapsis omega must be 120 degrees.
    # This gives a start true anomaly (at North Pole) of -30 deg
    # and an end true anomaly (at Equator) of 60 deg.
    nu_start_deg = -30.0
    nu_end_deg = 60.0

    # Convert degrees to radians for calculations
    nu_start_rad = math.radians(nu_start_deg)
    nu_end_rad = math.radians(nu_end_deg)

    # Step 4: Convert true anomalies to eccentric anomalies (E)
    # Formula: E = 2 * atan(sqrt((1-e)/(1+e)) * tan(nu/2))
    eccentricity_term = math.sqrt((1 - e) / (1 + e))
    E_start_rad = 2 * math.atan(eccentricity_term * math.tan(nu_start_rad / 2))
    E_end_rad = 2 * math.atan(eccentricity_term * math.tan(nu_end_rad / 2))

    # Step 5: Convert eccentric anomalies to mean anomalies (M)
    # Formula: M = E - e * sin(E)
    M_start_rad = E_start_rad - e * math.sin(E_start_rad)
    M_end_rad = E_end_rad - e * math.sin(E_end_rad)

    # Step 6: Calculate mean motion (n)
    # Formula: n = sqrt(mu / a^3)
    n = math.sqrt(mu_m3_s2 / a**3)

    # Step 7: Calculate the time of flight
    delta_M = M_end_rad - M_start_rad
    time_s = delta_M / n

    # Step 8: Format the output
    # Round the final time to the nearest 10 seconds
    time_rounded = int(round(time_s / 10) * 10)

    # Print the final equation with its numerical components
    print("The time of flight is calculated as t = (M_end - M_start) / n")
    print(f"t = ({M_end_rad:.6f} rad - ({M_start_rad:.6f} rad)) / {n:.8f} rad/s")
    print(f"t = {delta_M:.6f} rad / {n:.8f} rad/s")
    print(f"t = {time_s:.2f} seconds")
    print("\nRounding to the nearest 10 seconds:")
    print(f"Final Answer: {time_rounded}")
    
    # Final answer in the required format
    print(f"<<<{time_rounded}>>>")

solve_orbital_time()