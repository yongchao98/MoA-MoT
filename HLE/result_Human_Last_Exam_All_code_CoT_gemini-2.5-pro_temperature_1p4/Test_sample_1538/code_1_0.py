import math

def solve_orbit_time():
    """
    Calculates the time for a spacecraft to travel from the north pole
    through periapsis to the equatorial plane around Mercury.
    """
    # Step 1: Define constants and convert to SI units (meters, seconds)
    R_mercury_km = 2440
    h_p_km = 218
    h_a_km = 9982
    mu_mercury_m3_s2 = 2.2032e13

    R_mercury = R_mercury_km * 1000
    h_p = h_p_km * 1000
    h_a = h_a_km * 1000

    # Step 2: Calculate orbital elements
    r_p = R_mercury + h_p
    r_a = R_mercury + h_a

    a = (r_p + r_a) / 2
    e = (r_a - r_p) / (r_a + r_p)

    # Step 3: Determine true anomalies for the start and end points
    # As determined in the plan, the path corresponds to:
    nu1_deg = -30.0  # True anomaly at the North Pole (start point)
    nu2_deg = 60.0   # True anomaly at the Equator (end point)

    # Convert degrees to radians for calculations
    nu1_rad = math.radians(nu1_deg)
    nu2_rad = math.radians(nu2_deg)

    # Step 4: Calculate the time of flight
    
    # a) Calculate mean motion (n)
    n = math.sqrt(mu_mercury_m3_s2 / (a**3))

    # b) Calculate eccentric anomalies (E) from true anomalies (nu)
    # The relation is: tan(E/2) = sqrt((1-e)/(1+e)) * tan(nu/2)
    conv_factor = math.sqrt((1 - e) / (1 + e))
    E1_rad = 2 * math.atan(conv_factor * math.tan(nu1_rad / 2))
    E2_rad = 2 * math.atan(conv_factor * math.tan(nu2_rad / 2))

    # c) Calculate mean anomalies (M) from eccentric anomalies (E)
    # The relation is Kepler's Equation: M = E - e*sin(E)
    M1_rad = E1_rad - e * math.sin(E1_rad)
    M2_rad = E2_rad - e * math.sin(E2_rad)

    # d) Calculate the time of flight (delta_t)
    delta_t_seconds = (M2_rad - M1_rad) / n

    # Step 5: Round the final answer to the nearest 10 seconds
    final_time = int(round(delta_t_seconds / 10.0)) * 10
    
    # Print the detailed calculation steps and the final answer
    print("--- Orbital Parameters ---")
    print(f"Periapsis radius (r_p): {r_p:.0f} m")
    print(f"Apoapsis radius (r_a): {r_a:.0f} m")
    print(f"Semi-major axis (a): {a:.0f} m")
    print(f"Eccentricity (e): {e:.6f}")
    
    print("\n--- Time of Flight Calculation ---")
    print("Equation: t = (M_2 - M_1) / n")
    print(f"Mean Anomaly at start (M_1): {M1_rad:.6f} radians (for true anomaly {nu1_deg}°)")
    print(f"Mean Anomaly at end (M_2): {M2_rad:.6f} radians (for true anomaly {nu2_deg}°)")
    print(f"Mean motion (n): {n:.8f} rad/s")
    
    print("\nSubstituting the numbers into the final equation:")
    print(f"t = ({M2_rad:.6f} - ({M1_rad:.6f})) / {n:.8f}")
    print(f"Calculated time: {delta_t_seconds:.2f} seconds")

    print("\n--- Final Answer ---")
    print("The time it took for the spacecraft to travel is (rounded to the nearest 10 seconds):")
    print(final_time)

solve_orbit_time()