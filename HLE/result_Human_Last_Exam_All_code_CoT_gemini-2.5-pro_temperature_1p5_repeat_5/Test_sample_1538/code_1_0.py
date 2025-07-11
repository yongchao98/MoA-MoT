import math

def solve_orbit_time():
    """
    Calculates the time for a spacecraft to travel from the North Pole,
    through periapsis, to the equator of Mercury.
    """
    # Step 1: Given Information in SI units
    h_p = 218 * 1000  # Periapsis altitude in meters
    h_a = 9982 * 1000  # Apoapsis altitude in meters
    R_M = 2440 * 1000  # Radius of Mercury in meters
    mu = 2.2032e13  # Gravitational parameter of Mercury in m^3/s^2

    # Step 2: Calculate Orbital Parameters
    r_p = R_M + h_p
    r_a = R_M + h_a
    a = (r_p + r_a) / 2
    e = (r_a - r_p) / (r_a + r_p)

    # Step 3 & 4: Determine True Anomalies for the trajectory
    # The path from North Pole -> Periapsis -> Equator requires an argument of
    # periapsis (omega) of 120 degrees for prograde motion.
    # Start Point (North Pole, lat=90deg): nu_start = -30 deg
    # End Point (Equator, lat=0deg): nu_end = 60 deg
    nu_start_rad = math.radians(-30)
    nu_end_rad = math.radians(60)

    # Step 5: Calculate Time of Flight using Kepler's Equation

    # Function to get Eccentric Anomaly (E) from True Anomaly (nu)
    def get_eccentric_anomaly(nu_rad, e):
        term = math.sqrt((1 - e) / (1 + e)) * math.tan(nu_rad / 2)
        return 2 * math.atan(term)

    # Function to get Mean Anomaly (M) from Eccentric Anomaly (E)
    def get_mean_anomaly(E_rad, e):
        return E_rad - e * math.sin(E_rad)

    # Calculate E and M for start and end points
    E_start_rad = get_eccentric_anomaly(nu_start_rad, e)
    E_end_rad = get_eccentric_anomaly(nu_end_rad, e)

    M_start_rad = get_mean_anomaly(E_start_rad, e)
    M_end_rad = get_mean_anomaly(E_end_rad, e)

    # Calculate the mean motion term and the time difference
    mean_motion_term = math.sqrt(a**3 / mu)
    delta_M = M_end_rad - M_start_rad
    delta_t = mean_motion_term * delta_M

    # Print the final equation components
    print("The final equation for the time of flight is: delta_t = sqrt(a^3 / mu) * (M_end - M_start)")
    print(f"Value for sqrt(a^3 / mu): {mean_motion_term:.4f} s/rad")
    print(f"Value for a (semi-major axis): {a} m")
    print(f"Value for mu (gravitational parameter): {mu} m^3/s^2")
    print(f"Value for M_end: {M_end_rad:.4f} rad")
    print(f"Value for M_start: {M_start_rad:.4f} rad")
    print(f"Calculated time of flight (delta_t): {delta_t:.4f} s")
    
    # Step 6: Round to the nearest 10 seconds
    final_time = int(round(delta_t / 10) * 10)
    
    print(f"\nThe final answer rounded to the nearest 10 seconds is: {final_time}")


solve_orbit_time()