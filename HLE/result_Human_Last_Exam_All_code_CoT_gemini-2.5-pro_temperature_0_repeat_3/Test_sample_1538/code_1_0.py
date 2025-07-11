import math

def solve_orbit_time():
    """
    Calculates the time for a spacecraft to travel from the North Pole,
    through periapsis, to the equatorial plane of Mercury.
    """
    # Step 1: Define constants and convert to SI units
    R_M = 2440 * 1000  # Radius of Mercury in meters
    mu_M = 2.2032e13   # Standard gravitational parameter of Mercury in m^3/s^2
    h_p = 218 * 1000   # Periapsis altitude in meters
    h_a = 9982 * 1000  # Apoapsis altitude in meters

    print("--- Given Constants ---")
    print(f"Radius of Mercury (R): {R_M} m")
    print(f"Gravitational Parameter (μ): {mu_M} m^3/s^2")
    print(f"Periapsis Altitude (h_p): {h_p} m")
    print(f"Apoapsis Altitude (h_a): {h_a} m")
    print("-" * 25)

    # Step 2: Calculate orbital parameters
    r_p = R_M + h_p
    r_a = R_M + h_a
    a = (r_p + r_a) / 2
    e = (r_a - r_p) / (r_a + r_p)
    n = math.sqrt(mu_M / a**3)

    print("--- Calculated Orbital Parameters ---")
    print(f"Radius of Periapsis (r_p): {r_p} m")
    print(f"Radius of Apoapsis (r_a): {r_a} m")
    print(f"Semi-major Axis (a): {a} m")
    print(f"Eccentricity (e): {e}")
    print(f"Mean Motion (n): {n} rad/s")
    print("-" * 25)

    # Step 3: Determine true anomalies for start and end points
    # As determined by the geometry (polar orbit, periapsis at 60N lat):
    # Start point (North Pole) has true anomaly nu1 = -30 deg
    # End point (Equator) has true anomaly nu2 = 60 deg
    nu1_deg = -30.0
    nu2_deg = 60.0
    nu1_rad = math.radians(nu1_deg)
    nu2_rad = math.radians(nu2_deg)
    
    print("--- Flight Path Anomalies ---")
    print(f"Start True Anomaly (ν1): {nu1_deg} degrees")
    print(f"End True Anomaly (ν2): {nu2_deg} degrees")
    print("-" * 25)

    # Step 4: Function to calculate time from periapsis
    def get_time_from_periapsis(nu_rad, e, n):
        # Convert true anomaly (nu) to eccentric anomaly (E)
        E = 2 * math.atan2(math.sqrt(1 - e) * math.sin(nu_rad / 2), math.sqrt(1 + e) * math.cos(nu_rad / 2))
        # Calculate mean anomaly (M) from eccentric anomaly (E)
        M = E - e * math.sin(E)
        # Calculate time (t) from mean anomaly (M)
        t = M / n
        return t, E, M

    # Calculate times for start and end points
    t1, E1, M1 = get_time_from_periapsis(nu1_rad, e, n)
    t2, E2, M2 = get_time_from_periapsis(nu2_rad, e, n)

    print("--- Time of Flight Calculation ---")
    print(f"For Start Point (ν1 = {nu1_deg}°):")
    print(f"  Eccentric Anomaly (E1): {math.degrees(E1):.4f} degrees")
    print(f"  Mean Anomaly (M1): {math.degrees(M1):.4f} degrees")
    print(f"  Time from Periapsis (t1): {t1:.2f} s")
    
    print(f"\nFor End Point (ν2 = {nu2_deg}°):")
    print(f"  Eccentric Anomaly (E2): {math.degrees(E2):.4f} degrees")
    print(f"  Mean Anomaly (M2): {math.degrees(M2):.4f} degrees")
    print(f"  Time from Periapsis (t2): {t2:.2f} s")
    print("-" * 25)

    # Calculate total travel time
    delta_t = t2 - t1

    # Step 5: Round to the nearest 10 seconds
    final_time = int(round(delta_t / 10) * 10)

    print(f"Total travel time (Δt = t2 - t1): {delta_t:.2f} s")
    print(f"\nFinal Answer: The time it took for the spacecraft to travel, rounded to the nearest 10 seconds, is {final_time} seconds.")

    return final_time

if __name__ == '__main__':
    answer = solve_orbit_time()
    print(f"<<<{answer}>>>")
