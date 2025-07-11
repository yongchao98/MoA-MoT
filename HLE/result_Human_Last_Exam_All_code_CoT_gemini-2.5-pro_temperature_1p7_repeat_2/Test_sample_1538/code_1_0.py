import math

def solve_orbital_time():
    """
    Calculates the time for a spacecraft to travel from the North Pole,
    through periapsis, to the equatorial plane around Mercury.
    """
    # Step 1: Define constants and calculate basic orbital parameters
    h_p = 218 * 1000  # Periapsis altitude in meters
    h_a = 9982 * 1000  # Apoapsis altitude in meters
    R_M = 2440 * 1000  # Radius of Mercury in meters
    mu = 2.2032e13  # Standard gravitational parameter of Mercury in m^3/s^2

    # Calculate periapsis and apoapsis radii from the center of Mercury
    r_p = R_M + h_p
    r_a = R_M + h_a

    # Calculate semi-major axis (a) and eccentricity (e)
    a = (r_p + r_a) / 2
    e = (r_a - r_p) / (r_a + r_p)

    # Step 2: Determine true anomalies for start and end points
    # For a polar orbit (i=90 deg), latitude relates to true anomaly (nu)
    # and argument of periapsis (omega) by: sin(latitude) = sin(omega + nu).
    # At periapsis (nu=0), latitude=60N. sin(60) = sin(omega).
    # The path is "through periapsis" from a higher to a lower latitude, which
    # implies omega must be 120 degrees for the geometry to work.
    # Start (North Pole, 90N): sin(90) = 1 = sin(120 + nu1) -> 120 + nu1 = 90 -> nu1 = -30 deg.
    # End (Equator, 0): sin(0) = 0 = sin(120 + nu2) -> 120 + nu2 = 180 -> nu2 = 60 deg.
    nu1_deg = -30.0
    nu2_deg = 60.0
    nu1_rad = math.radians(nu1_deg)
    nu2_rad = math.radians(nu2_deg)

    # Step 3: Convert true anomalies to eccentric anomalies (E)
    def true_to_eccentric(nu_rad, e):
        """Converts true anomaly (radians) to eccentric anomaly (radians)."""
        tan_E_over_2 = math.sqrt((1 - e) / (1 + e)) * math.tan(nu_rad / 2)
        return 2 * math.atan(tan_E_over_2)

    E1_rad = true_to_eccentric(nu1_rad, e)
    E2_rad = true_to_eccentric(nu2_rad, e)

    # Step 4: Use Kepler's Equation to find the time of flight
    # Calculate mean anomalies (M) from eccentric anomalies (E)
    M1_rad = E1_rad - e * math.sin(E1_rad)
    M2_rad = E2_rad - e * math.sin(E2_rad)

    # Calculate mean motion (n)
    n = math.sqrt(mu / (a**3))

    # Calculate time of flight (delta_t)
    delta_t_seconds = (M2_rad - M1_rad) / n

    # Step 5: Round the final answer to the nearest 10 seconds
    final_answer = int(round(delta_t_seconds, -1))
    
    # Print the values used in the final equation
    print("--- Time of Flight Calculation Breakdown ---")
    print(f"The final time is calculated using the equation: Δt = (M2 - M1) / n")
    print("\nThe values for each number in this equation are:")
    print(f"a (semi-major axis) = {a} m")
    print(f"μ (gravitational parameter) = {mu} m^3/s^2")
    print(f"M1 (start mean anomaly) = {M1_rad} rad")
    print(f"M2 (end mean anomaly) = {M2_rad} rad")
    print(f"n (mean motion) = {n} rad/s")
    
    print("\nPlugging in the numbers:")
    print(f"Δt = ({M2_rad} - ({M1_rad})) / {n}")
    print(f"Δt = {M2_rad - M1_rad} / {n}")

    print(f"\nUnrounded Time of Flight = {delta_t_seconds} seconds")
    print(f"Final Answer (rounded to the nearest 10s) = {final_answer}")

solve_orbital_time()