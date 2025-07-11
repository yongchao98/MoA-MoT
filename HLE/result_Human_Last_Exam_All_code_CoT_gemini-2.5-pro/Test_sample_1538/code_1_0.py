import math

def solve_orbital_time():
    """
    Calculates the time it took for a spacecraft to travel from the north pole,
    through periapsis, to the equatorial plane of Mercury.
    """
    # Step 1: Define constants and convert units to meters
    R_M = 2440 * 1000  # Radius of Mercury in meters
    mu = 2.2032e13     # Standard gravitational parameter of Mercury in m^3/s^2
    h_p = 218 * 1000   # Periapsis altitude in meters
    h_a = 9982 * 1000  # Apoapsis altitude in meters

    # Step 2: Calculate orbital parameters
    r_p = R_M + h_p  # Radius of periapsis
    r_a = R_M + h_a  # Radius of apoapsis
    a = (r_p + r_a) / 2  # Semi-major axis
    e = (r_a - r_p) / (r_a + r_p)  # Eccentricity

    # Step 3: Determine true anomalies for the start and end points
    # For a polar orbit (i=90 deg), sin(latitude) = sin(argument_of_periapsis + true_anomaly)
    # At periapsis (nu=0), latitude is 60 deg N. So, sin(60) = sin(omega). Thus, omega = 60 deg.
    omega_rad = math.radians(60)

    # Start point: North Pole (latitude = 90 deg N), before periapsis
    # sin(90) = sin(60 + nu_np) => 1 = sin(60 + nu_np) => 60 + nu_np = 90 => nu_np = 30 deg
    # Since the point is before periapsis, the true anomaly is negative.
    nu_np_deg = -30.0

    # End point: Equatorial Plane (latitude = 0 deg), after periapsis
    # sin(0) = sin(60 + nu_eq) => 0 = sin(60 + nu_eq) => 60 + nu_eq = 180 (first crossing after periapsis)
    # => nu_eq = 120 deg
    nu_eq_deg = 120.0

    # Step 4: Create a function to calculate time from periapsis passage
    def get_time_from_periapsis(nu_deg, a_val, e_val, mu_val):
        """Calculates time from periapsis for a given true anomaly."""
        nu_rad = math.radians(nu_deg)

        # Convert true anomaly to eccentric anomaly (E)
        # Handle atan2 for correct quadrant
        E_rad = 2 * math.atan2(math.sqrt(1 - e_val) * math.sin(nu_rad / 2),
                               math.sqrt(1 + e_val) * math.cos(nu_rad / 2))

        # Calculate mean anomaly (M) using Kepler's Equation
        M_rad = E_rad - e_val * math.sin(E_rad)

        # Calculate mean motion (n)
        n = math.sqrt(mu_val / a_val**3)

        # Calculate time from periapsis
        t = M_rad / n
        return t

    # Calculate time for start and end points
    t_np = get_time_from_periapsis(nu_np_deg, a, e, mu)
    t_eq = get_time_from_periapsis(nu_eq_deg, a, e, mu)

    # The total time of flight is the difference
    time_of_flight = t_eq - t_np

    # Step 5: Round the final answer to the nearest 10 seconds
    final_time = int(round(time_of_flight / 10.0) * 10)
    
    print(f"The time it took for the spacecraft to travel is {final_time} seconds.")

solve_orbital_time()
<<<3560>>>