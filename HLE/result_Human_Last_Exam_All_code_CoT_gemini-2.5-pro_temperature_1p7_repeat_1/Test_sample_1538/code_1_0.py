import math

def calculate_time_of_flight():
    """
    Calculates the time it took for the spacecraft to travel from the north pole,
    through periapsis, to the equatorial plane of Mercury.
    """
    # Step 1: Define constants and calculate orbital parameters
    R_mercury = 2440 * 1000  # Radius of Mercury in meters
    mu_mercury = 2.2032e13    # Gravitational parameter of Mercury in m^3/s^2
    h_p = 218 * 1000         # Periapsis altitude in meters
    h_a = 9982 * 1000        # Apoapsis altitude in meters

    r_p = R_mercury + h_p  # Periapsis radius
    r_a = R_mercury + h_a  # Apoapsis radius

    a = (r_p + r_a) / 2.0  # Semi-major axis
    e = (r_a - r_p) / (r_a + r_p)  # Eccentricity
    n = math.sqrt(mu_mercury / a**3)  # Mean motion in rad/s

    # Step 2: Define start and end true anomalies based on the problem statement
    # As derived from the orbit geometry (polar orbit, periapsis at 60N):
    # Start at North Pole (90N) corresponds to nu_1 = -30 degrees.
    # End at Equator (0N) corresponds to nu_2 = 60 degrees.
    v1_rad = math.radians(-30)
    v2_rad = math.radians(60)

    # Step 3: Calculate Mean Anomaly (M) for both points using a helper function
    def get_mean_anomaly(v_rad, ecc):
        """Converts true anomaly (v) to mean anomaly (M)."""
        # Convert true anomaly to eccentric anomaly (E)
        tan_E_half = math.sqrt((1 - ecc) / (1 + ecc)) * math.tan(v_rad / 2.0)
        E = 2 * math.atan(tan_E_half)
        # Convert eccentric anomaly to mean anomaly (M) via Kepler's Equation
        M = E - ecc * math.sin(E)
        return M

    M1 = get_mean_anomaly(v1_rad, e)
    M2 = get_mean_anomaly(v2_rad, e)

    # Step 4: Calculate the time of flight
    time_of_flight = (M2 - M1) / n

    # Step 5: Round to the nearest 10 seconds
    final_time = int(round(time_of_flight / 10.0) * 10)
    
    # Output the final equation with the calculated numbers
    print(f"The time of flight is calculated from the equation: Time = (M2 - M1) / n")
    print(f"Time = ({M2} - ({M1})) / {n}")

    # Output the final rounded answer
    print(f"The calculated time of flight is {final_time} seconds.")
    
    # Return the answer in the required format
    print(f"<<<{final_time}>>>")

# Execute the calculation
calculate_time_of_flight()