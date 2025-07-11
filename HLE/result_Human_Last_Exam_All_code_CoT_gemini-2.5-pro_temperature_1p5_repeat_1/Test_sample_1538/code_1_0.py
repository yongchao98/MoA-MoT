import math

def solve_orbital_time():
    """
    Calculates the time it took for a spacecraft to travel from the north pole
    through periapsis to the equatorial plane of Mercury.
    """

    # 1. Define Constants and Convert Units
    R_mercury = 2440 * 1000  # Radius of Mercury in meters
    mu_mercury = 2.2032e13  # Standard gravitational parameter of Mercury in m^3/s^2
    h_p = 218 * 1000  # Periapsis altitude in meters
    h_a = 9982 * 1000 # Apoapsis altitude in meters

    # 2. Calculate Orbital Radii
    r_p = R_mercury + h_p  # Radius of periapsis
    r_a = R_mercury + h_a  # Radius of apoapsis

    # 3. Determine Key Orbital Parameters
    a = (r_p + r_a) / 2  # Semi-major axis
    e = (r_a - r_p) / (r_a + r_p)  # Eccentricity
    n = math.sqrt(mu_mercury / a**3)  # Mean motion

    # 4. Determine True Anomalies
    # The problem states the spacecraft travels from the north pole (90 N) to the equator (0).
    # Argument of periapsis (omega) is 60 deg, since periapsis is at 60 N.
    # For a polar orbit, latitude = omega + true_anomaly.
    nu1_deg = 90.0 - 60.0  # True anomaly at start point (North Pole)
    nu2_deg = 0.0 - 60.0   # True anomaly at end point (Equator)

    nu1_rad = math.radians(nu1_deg)
    nu2_rad = math.radians(nu2_deg)

    # 5. Function to Calculate Mean Anomaly from True Anomaly
    def get_mean_anomaly(nu_rad, eccentricity):
        """Calculates mean anomaly from true anomaly."""
        # Calculate eccentric anomaly (E) from true anomaly (nu)
        tan_E_over_2 = math.sqrt((1 - eccentricity) / (1 + eccentricity)) * math.tan(nu_rad / 2)
        E_rad = 2 * math.atan(tan_E_over_2)
        # Calculate mean anomaly (M) from eccentric anomaly (E) using Kepler's Equation
        M_rad = E_rad - eccentricity * math.sin(E_rad)
        return M_rad

    # Calculate mean anomalies for start and end points
    M1_rad = get_mean_anomaly(nu1_rad, e)
    M2_rad = get_mean_anomaly(nu2_rad, e)

    # 6. Calculate Time of Flight (delta_t)
    delta_t = (M1_rad - M2_rad) / n

    # 7. Finalize and Print
    rounded_delta_t = int(round(delta_t / 10.0) * 10)

    print("--- Orbital Parameters ---")
    print(f"Periapsis Radius (r_p): {r_p} m")
    print(f"Apoapsis Radius (r_a): {r_a} m")
    print(f"Semi-major Axis (a): {a} m")
    print(f"Eccentricity (e): {e}")
    print(f"Mean Motion (n): {n} rad/s")
    print("\n--- Time of Flight Calculation ---")
    print(f"Start True Anomaly (nu1): {nu1_deg} degrees")
    print(f"End True Anomaly (nu2): {nu2_deg} degrees")
    print(f"Start Mean Anomaly (M1): {M1_rad} rad")
    print(f"End Mean Anomaly (M2): {M2_rad} rad")
    print("\n--- Final Equation ---")
    print(f"Time (s) = (M1 - M2) / n")
    print(f"Time (s) = ({M1_rad} - ({M2_rad})) / {n}")
    print(f"Calculated Time of Flight: {delta_t} s")
    print(f"\nRounded Answer (to nearest 10 seconds): {rounded_delta_t}")


solve_orbital_time()
<<<1270>>>