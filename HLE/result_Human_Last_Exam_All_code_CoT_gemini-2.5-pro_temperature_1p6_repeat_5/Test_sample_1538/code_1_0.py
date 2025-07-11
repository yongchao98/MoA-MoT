import math

def solve_orbital_time():
    """
    Calculates the time of flight for a spacecraft in a polar orbit around Mercury.
    """
    # Step 1: Define constants and orbital parameters in SI units
    R_mercury = 2440 * 1000  # Radius of Mercury in meters
    mu_mercury = 2.2032 * 10**13  # Gravitational parameter of Mercury in m^3/s^2
    h_p = 218 * 1000  # Periapsis altitude in meters
    h_a = 9982 * 1000  # Apoapsis altitude in meters

    # Step 2: Calculate primary orbital characteristics
    r_p = R_mercury + h_p  # Radius of periapsis
    r_a = R_mercury + h_a  # Radius of apoapsis
    a = (r_p + r_a) / 2  # Semi-major axis
    e = (r_a - r_p) / (r_a + r_p)  # Eccentricity

    # Step 3: Determine the true anomalies for the start and end points
    # The trajectory goes from the North Pole (u=90 deg) through periapsis to the Equator (u=180 deg).
    # Periapsis is at 60 deg latitude. For a polar orbit, sin(lat) = sin(u).
    # For periapsis to be on this path, its argument of latitude must be 120 deg.
    # At periapsis, true anomaly (theta) is 0, so the argument of periapsis (omega) is 120 deg.
    omega_rad = math.radians(120)

    # Start point: North Pole (u = 90 deg)
    u1_rad = math.radians(90)
    theta1_rad = u1_rad - omega_rad  # True anomaly at the start

    # End point: Equatorial plane (u = 180 deg)
    u2_rad = math.radians(180)
    theta2_rad = u2_rad - omega_rad  # True anomaly at the end

    # Step 4: Convert true anomalies (theta) to eccentric anomalies (E)
    def true_to_eccentric_anomaly(theta, eccentricity):
        """Converts true anomaly to eccentric anomaly."""
        return 2 * math.atan(math.sqrt((1 - eccentricity) / (1 + eccentricity)) * math.tan(theta / 2))

    E1_rad = true_to_eccentric_anomaly(theta1_rad, e)
    E2_rad = true_to_eccentric_anomaly(theta2_rad, e)

    # Step 5: Calculate the corresponding mean anomalies (M)
    M1_rad = E1_rad - e * math.sin(E1_rad)
    M2_rad = E2_rad - e * math.sin(E2_rad)

    # Step 6: Calculate the mean motion (n)
    n = math.sqrt(mu_mercury / (a**3))

    # Step 7: Calculate the time of flight (delta_t)
    delta_t_sec = (M2_rad - M1_rad) / n

    # Output the numbers used in the final equation as requested
    print(f"To find the time of flight (Δt), we use the equation: Δt = (M₂ - M₁) / n")
    print("\n--- Values in Final Equation ---")
    print(f"Mean Anomaly at start (M₁): {M1_rad} rad")
    print(f"Mean Anomaly at end (M₂): {M2_rad} rad")
    print(f"Mean Motion (n): {n} rad/s")
    print("---")
    
    print(f"\nFinal calculation is: Δt = ({M2_rad} - ({M1_rad})) / {n}")
    print(f"Result (unrounded): {delta_t_sec} seconds")


    # Step 8: Round the final answer to the nearest 10 seconds
    final_time = int(round(delta_t_sec, -1))
    
    print(f"\nRounded to the nearest 10 seconds, the final answer is:")
    print(final_time)
    
    return final_time

if __name__ == '__main__':
    solve_orbital_time()
