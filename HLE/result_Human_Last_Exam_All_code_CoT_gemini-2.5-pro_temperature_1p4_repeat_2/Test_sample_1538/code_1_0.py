import math

def solve_flight_time():
    """
    Calculates the time it took for the spacecraft to travel from the north pole 
    through periapsis to the equatorial plane of Mercury.
    """
    # Step 1: Define constants in SI units (meters, seconds).
    R_M = 2440 * 1000   # Radius of Mercury in meters
    mu = 2.2032e13      # Standard gravitational parameter of Mercury in m^3/s^2
    h_p = 218 * 1000    # Periapsis altitude in meters
    h_a = 9982 * 1000   # Apoapsis altitude in meters

    # Step 2: Calculate radii and orbital elements.
    r_p = R_M + h_p
    r_a = R_M + h_a
    a = (r_p + r_a) / 2
    e = (r_a - r_p) / (r_a + r_p)

    # Step 3: Define the start and end true anomalies in radians.
    # As determined in the plan, the path is from theta = -30 deg to theta = 60 deg.
    theta_1_rad = math.radians(-30.0)
    theta_2_rad = math.radians(60.0)

    # Step 4: Calculate time of flight.
    
    # Helper function to calculate mean anomaly (M) from true anomaly (theta)
    def get_mean_anomaly(true_anomaly_rad, eccentricity):
        """Calculates mean anomaly from true anomaly via eccentric anomaly."""
        # Calculate eccentric anomaly (E) from true anomaly (theta)
        # Using E = atan2(sqrt(1-e^2)*sin(theta), e + cos(theta)) for robustness
        E = math.atan2(math.sqrt(1 - eccentricity**2) * math.sin(true_anomaly_rad),
                       eccentricity + math.cos(true_anomaly_rad))
        
        # Calculate mean anomaly (M) using Kepler's Equation: M = E - e*sin(E)
        M = E - eccentricity * math.sin(E)
        return M

    # Calculate mean anomalies for the start and end points.
    M1 = get_mean_anomaly(theta_1_rad, e)
    M2 = get_mean_anomaly(theta_2_rad, e)

    # Calculate mean motion (n).
    n = math.sqrt(mu / a**3)
    
    # The final equation for time of flight is: time = (M2 - M1) / n
    # We now output each number used in this final equation.
    print(f"M2 (End Mean Anomaly in radians): {M2}")
    print(f"M1 (Start Mean Anomaly in radians): {M1}")
    print(f"n (Mean Motion in rad/s): {n}")

    # Calculate the time of flight in seconds.
    time_delta_seconds = (M2 - M1) / n

    # Step 5: Round the final answer to the nearest 10 seconds.
    rounded_time = int(round(time_delta_seconds / 10.0) * 10)

    print(f"\nFinal calculated time in seconds (rounded to nearest 10):")
    print(rounded_time)


solve_flight_time()