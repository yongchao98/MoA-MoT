import math

def solve_orbit_time():
    """
    Calculates the time for a spacecraft to travel from the North Pole,
    through periapsis, to the equatorial plane around Mercury.
    """
    # Step 1: Define given constants and calculate orbital parameters.
    # Altitudes and radius are converted to meters.
    h_p = 218 * 1000  # Periapsis altitude in meters
    h_a = 9982 * 1000 # Apoapsis altitude in meters
    R_M = 2440 * 1000 # Radius of Mercury in meters
    mu = 2.2032e13    # Gravitational parameter of Mercury in m^3/s^2

    # Calculate radii from the center of Mercury
    r_p = R_M + h_p
    r_a = R_M + h_a

    # Calculate semi-major axis (a) and eccentricity (e)
    a = (r_p + r_a) / 2
    e = (r_a - r_p) / (r_a + r_p)

    # Calculate mean motion (n) in radians per second
    n = math.sqrt(mu / a**3)

    # Step 2: Determine the true anomalies for the start and end points.
    # The relationship for a polar orbit (inclination i=90 deg) is:
    # sin(latitude) = sin(argument_of_periapsis + true_anomaly)
    # sin(lat) = sin(omega + nu)
    # At periapsis, nu = 0 and latitude = 60 deg.
    # sin(60) = sin(omega). This gives omega = 60 or 120 deg.
    # To have the North Pole passage *before* periapsis (nu < 0), we must choose omega = 120 deg.
    # Start point: North Pole (latitude = 90 deg)
    # sin(90) = sin(120 + nu1) => 1 = sin(120 + nu1) => 120 + nu1 = 90 => nu1 = -30 deg.
    # End point: Equatorial Plane (latitude = 0 deg) after periapsis
    # sin(0) = sin(120 + nu2) => 0 = sin(120 + nu2) => 120 + nu2 = 180 => nu2 = 60 deg.
    
    nu1_rad = math.radians(-30.0)
    nu2_rad = math.radians(60.0)

    # Step 3: Convert true anomalies to eccentric anomalies (E).
    # tan(E/2) = sqrt((1-e)/(1+e)) * tan(nu/2)
    def true_to_eccentric(nu_rad, eccentricity):
        tan_E_half = math.sqrt((1 - eccentricity) / (1 + eccentricity)) * math.tan(nu_rad / 2)
        return 2 * math.atan(tan_E_half)

    E1_rad = true_to_eccentric(nu1_rad, e)
    E2_rad = true_to_eccentric(nu2_rad, e)

    # Step 4: Convert eccentric anomalies to mean anomalies (M) using Kepler's Equation.
    # M = E - e*sin(E)
    def eccentric_to_mean(E_rad, eccentricity):
        return E_rad - eccentricity * math.sin(E_rad)

    M1_rad = eccentric_to_mean(E1_rad, e)
    M2_rad = eccentric_to_mean(E2_rad, e)

    # Step 5: Calculate the time of flight (delta_t).
    # delta_t = (M2 - M1) / n
    delta_t_seconds = (M2_rad - M1_rad) / n

    # Step 6: Round the final answer to the nearest 10 seconds.
    final_time = round(delta_t_seconds / 10) * 10
    
    # Print the details of the final calculation
    print(f"The final calculation for the time of flight is (M2 - M1) / n:")
    print(f"Time = ({M2_rad:.6f} rad - ({M1_rad:.6f} rad)) / {n:.8f} rad/s")
    print(f"Time = {delta_t_seconds:.2f} s")
    print(f"\nRounding to the nearest 10 seconds, the final answer is: {final_time}")
    
    return final_time

# Run the solver and print the final answer in the required format.
final_answer = solve_orbit_time()
print(f"<<<{final_answer}>>>")
