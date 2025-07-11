import numpy as np

def solve_orbital_time():
    """
    Calculates the time for a spacecraft to travel from the north pole
    through periapsis to the equator of Mercury.
    """
    # Step 1: Define constants and convert to SI units
    R_mercury = 2440 * 1000       # Radius of Mercury in meters
    mu_mercury = 2.2032e13      # Gravitational parameter of Mercury in m^3/s^2
    h_p = 218 * 1000            # Periapsis altitude in meters
    h_a = 9982 * 1000           # Apoapsis altitude in meters

    # Step 2: Calculate orbital parameters
    r_p = R_mercury + h_p
    r_a = R_mercury + h_a
    a = (r_p + r_a) / 2
    e = (r_a - r_p) / (r_a + r_p)
    n = np.sqrt(mu_mercury / a**3) # Mean motion in rad/s

    # Step 3: Determine true anomalies (in radians)
    # Start point: North Pole (lambda=90 deg). nu_1 = 90 - 60 = 30 deg
    nu1_rad = np.radians(30.0)
    # End point: Equator (lambda=0 deg). nu_2 = 0 - 60 = -60 deg
    nu2_rad = np.radians(-60.0)

    # Step 4: Calculate eccentric and mean anomalies
    # For start point (nu1)
    tan_E1_over_2 = np.sqrt((1 - e) / (1 + e)) * np.tan(nu1_rad / 2)
    E1_rad = 2 * np.arctan(tan_E1_over_2)
    M1_rad = E1_rad - e * np.sin(E1_rad)

    # For end point (nu2)
    tan_E2_over_2 = np.sqrt((1 - e) / (1 + e)) * np.tan(nu2_rad / 2)
    E2_rad = 2 * np.arctan(tan_E2_over_2)
    M2_rad = E2_rad - e * np.sin(E2_rad)

    # Step 5: Calculate time of flight
    # The time taken is t = (M1 - M2) / n
    # where M1 corresponds to the time from periapsis to nu=30deg
    # and M2 corresponds to the time from periapsis to nu=-60deg (which is negative)
    # The total time from nu=30, through nu=0, to nu=-60 is the sum of magnitudes of these times.
    delta_t = (M1_rad - M2_rad) / n

    # Print the numbers used in the final equation
    print(f"The total time is calculated using the formula: (M1 - M2) / n")
    print(f"M1 (Mean Anomaly at start point): {M1_rad:.6f} radians")
    print(f"M2 (Mean Anomaly at end point): {M2_rad:.6f} radians")
    print(f"n (Mean Motion): {n:.6f} rad/s")
    print(f"Time in seconds = ({M1_rad:.6f} - ({M2_rad:.6f})) / {n:.6f}")
    
    # Round the final answer to the nearest 10 seconds
    final_time = int(round(delta_t, -1))

    print(f"\nThe calculated time is {delta_t:.2f} seconds.")
    print(f"The final answer rounded to the nearest 10 seconds is: {final_time}")
    
    # Return final answer in specified format
    return final_time

if __name__ == '__main__':
    answer = solve_orbital_time()
    print(f"\n<<< {answer} >>>")
