import math

def solve_flight_time():
    """
    Calculates the time of flight for a spacecraft in a specific orbital segment around Mercury.
    """
    # Given constants
    R_M = 2440 * 1000   # Radius of Mercury in meters
    h_p = 218 * 1000    # Periapsis altitude in meters
    h_a = 9982 * 1000   # Apoapsis altitude in meters
    mu = 2.2032e13      # Standard gravitational parameter of Mercury in m^3/s^2

    # Step 1: Calculate orbital parameters from the given altitudes and constants.
    r_p = R_M + h_p
    r_a = R_M + h_a
    a = (r_p + r_a) / 2
    e = (r_a - r_p) / (r_a + r_p)
    n = math.sqrt(mu / a**3) # Mean motion

    # Step 2: Define true anomalies for the start and end points in degrees.
    # As derived in the plan, start is at the North Pole (nu=330 deg) and
    # end is at the Equator after periapsis (nu=60 deg).
    nu1_deg = 330.0
    nu2_deg = 60.0

    # Convert degrees to radians for calculations
    nu1_rad = math.radians(nu1_deg)
    nu2_rad = math.radians(nu2_deg)

    # Step 3: Define helper functions to calculate eccentric and mean anomalies.
    def get_eccentric_anomaly(nu_rad, eccentricity):
        """Calculates eccentric anomaly (E) from true anomaly (nu)."""
        E_rad = math.atan2(math.sqrt(1 - eccentricity**2) * math.sin(nu_rad), eccentricity + math.cos(nu_rad))
        # atan2 returns in (-pi, pi], adjust to [0, 2*pi) for consistency
        if E_rad < 0:
            E_rad += 2 * math.pi
        return E_rad

    def get_mean_anomaly(E_rad, eccentricity):
        """Calculates mean anomaly (M) from eccentric anomaly (E)."""
        return E_rad - eccentricity * math.sin(E_rad)

    # Step 4: Calculate the eccentric and mean anomalies for the start and end points.
    E1_rad = get_eccentric_anomaly(nu1_rad, e)
    M1_rad = get_mean_anomaly(E1_rad, e)

    E2_rad = get_eccentric_anomaly(nu2_rad, e)
    M2_rad = get_mean_anomaly(E2_rad, e)

    # Step 5: Calculate the total time of flight.
    # The formula handles the flight path crossing periapsis.
    delta_t = (2 * math.pi + M2_rad - M1_rad) / n

    # Step 6: Print the components of the final calculation and the result.
    print("The time of flight (delta_t) is calculated using the formula: delta_t = (2*pi + M_2 - M_1) / n")
    print("\nSubstituting the calculated values:")
    
    term_2pi = 2 * math.pi
    
    print(f"delta_t = ({term_2pi:.6f} + {M2_rad:.6f} - {M1_rad:.6f}) / {n:.8f}")
    
    numerator = term_2pi + M2_rad - M1_rad
    print(f"delta_t = {numerator:.6f} / {n:.8f}")
    
    # Final result in seconds
    print(f"\nCalculated time of flight: {delta_t:.2f} seconds")
    
    # Round to the nearest 10 seconds as required.
    final_answer = int(round(delta_t, -1))
    
    print(f"The final answer rounded to the nearest 10 seconds is: {final_answer}")
    
    # Output the final answer in the required format
    print(f"\n<<<{final_answer}>>>")

solve_flight_time()