import math

def solve_orbital_time():
    """
    Calculates the time for a spacecraft to travel from Mercury's north pole,
    through periapsis, to the equatorial plane.
    """
    # Step 1: Define constants and initial parameters in SI units
    R_M = 2440 * 1000  # Radius of Mercury in meters
    h_p = 218 * 1000   # Periapsis altitude in meters
    h_a = 9982 * 1000  # Apoapsis altitude in meters
    mu = 2.2032e13     # Gravitational parameter of Mercury in m^3/s^2

    # Step 2: Calculate orbital parameters
    r_p = R_M + h_p
    r_a = R_M + h_a
    a = (r_p + r_a) / 2
    e = (r_a - r_p) / (r_a + r_p)

    # Step 3: Determine true anomalies for start and end points
    # From the problem geometry (polar orbit, periapsis at 60N, path is
    # N.Pole -> Periapsis -> Equator), we deduce:
    # True anomaly at the start point (North Pole): nu_1 = -30 degrees
    # True anomaly at the end point (Equator): nu_2 = 60 degrees
    nu1_rad = math.radians(-30)
    nu2_rad = math.radians(60)

    # Step 4: Solve Kepler's problem for time of flight
    # Function to calculate eccentric anomaly E from true anomaly nu
    def get_E(nu, e):
        return math.atan2(math.sqrt(1 - e**2) * math.sin(nu), e + math.cos(nu))

    # Function to calculate mean anomaly M from eccentric anomaly E
    def get_M(E, e):
        return E - e * math.sin(E)

    # Calculate eccentric and mean anomalies for the start point
    E1 = get_E(nu1_rad, e)
    M1 = get_M(E1, e)

    # Calculate eccentric and mean anomalies for the end point
    E2 = get_E(nu2_rad, e)
    M2 = get_M(E2, e)

    # Calculate mean motion (n)
    n = math.sqrt(mu / a**3)

    # Calculate the total time of flight
    delta_t = (M2 - M1) / n
    
    # Step 5: Round the result to the nearest 10 seconds
    final_time = int(round(delta_t / 10.0) * 10)

    # The problem asks to output the numbers in the final equation.
    # The final calculation is Time = (Mean Anomaly 2 - Mean Anomaly 1) / Mean Motion
    print("This problem is solved by finding the time of flight between two points in an orbit.")
    print("The final calculation is Time = (M2 - M1) / n")
    print(f"M1 (Mean Anomaly at North Pole) = {M1:.6f} radians")
    print(f"M2 (Mean Anomaly at Equator) = {M2:.6f} radians")
    print(f"n (Mean Motion) = {n:.6f} rad/s")
    print(f"Time = ({M2:.6f} - ({M1:.6f})) / {n:.6f}")
    print(f"Calculated time in seconds = {delta_t:.2f}")
    print(f"Time rounded to the nearest 10 seconds = {final_time}")
    
    # The final answer format is specified at the end.
    print("<<<1270>>>")

solve_orbital_time()