import math

def solve_orbital_time():
    """
    Calculates the time for a spacecraft to travel from the North Pole,
    through periapsis, to the equatorial plane of Mercury.
    """
    # Step 1: Define constants and convert units
    R_M = 2440 * 1000  # Radius of Mercury in meters
    h_p = 218 * 1000   # Periapsis altitude in meters
    h_a = 9982 * 1000  # Apoapsis altitude in meters
    mu = 2.2032e13     # Standard gravitational parameter of Mercury in m^3/s^2

    print("--- Given Parameters ---")
    print(f"Radius of Mercury (R_M): {R_M / 1000:.0f} km")
    print(f"Periapsis altitude (h_p): {h_p / 1000:.0f} km")
    print(f"Apoapsis altitude (h_a): {h_a / 1000:.0f} km")
    print(f"Standard gravitational parameter (mu): {mu:.4e} m^3 s^-2")
    print("-" * 26)

    # Step 2: Calculate orbital parameters
    r_p = R_M + h_p  # Radius of periapsis
    r_a = R_M + h_a  # Radius of apoapsis
    a = (r_p + r_a) / 2  # Semi-major axis
    e = (r_a - r_p) / (r_a + r_p)  # Eccentricity

    print("--- Calculated Orbital Parameters ---")
    print(f"Semi-major axis (a): {a:.0f} m")
    print(f"Eccentricity (e): {e:.6f}")
    print("-" * 35)

    # Step 3: Determine true anomalies
    # The path is North Pole (phi=90) -> Periapsis (nu=0) -> Equator (phi=0).
    # This requires nu_NP < 0 and nu_EQ > 0.
    # From sin(phi) = sin(omega + nu), we determine omega = 120 degrees.
    # Start (North Pole): sin(90) = sin(120 + nu1) -> 120 + nu1 = 90 -> nu1 = -30 deg
    # End (Equator): sin(0) = sin(120 + nu2) -> 120 + nu2 = 180 -> nu2 = 60 deg
    nu1_deg = -30.0
    nu2_deg = 60.0

    print("--- True Anomaly Calculation ---")
    print(f"Start true anomaly (nu1): {nu1_deg:.1f} degrees")
    print(f"End true anomaly (nu2): {nu2_deg:.1f} degrees")
    print("-" * 32)
    
    # Step 4: Calculate time of flight
    def get_mean_anomaly(nu_deg, e):
        """Calculates mean anomaly from true anomaly."""
        nu_rad = math.radians(nu_deg)
        # Convert true anomaly to eccentric anomaly
        tan_E_half = math.sqrt((1 - e) / (1 + e)) * math.tan(nu_rad / 2)
        E = 2 * math.atan(tan_E_half)
        # Convert eccentric anomaly to mean anomaly (Kepler's Equation)
        M = E - e * math.sin(E)
        return M

    M1 = get_mean_anomaly(nu1_deg, e)
    M2 = get_mean_anomaly(nu2_deg, e)
    
    delta_M = M2 - M1
    
    # Calculate time of flight
    time_of_flight = math.sqrt(a**3 / mu) * delta_M

    print("--- Time of Flight Calculation ---")
    print("Equation: delta_t = (M2 - M1) * sqrt(a^3 / mu)")
    print(f"M1 (radians): {M1:.6f}")
    print(f"M2 (radians): {M2:.6f}")
    print(f"delta_M (radians): {delta_M:.6f}")
    print(f"Term sqrt(a^3/mu) (s): {math.sqrt(a**3 / mu):.2f}")
    
    print("\n--- Final Equation with Values ---")
    print(f"Time (s) = ({M2:.6f} - ({M1:.6f})) * sqrt({a:.0f}^3 / {mu:.4e})")
    print(f"Time (s) = {time_of_flight:.2f}")

    # Round to the nearest 10 seconds
    final_answer = int(round(time_of_flight, -1))

    print("\n--- Final Answer ---")
    print(f"The time it took for the spacecraft to travel is {final_answer} seconds (rounded to the nearest 10).")

solve_orbital_time()