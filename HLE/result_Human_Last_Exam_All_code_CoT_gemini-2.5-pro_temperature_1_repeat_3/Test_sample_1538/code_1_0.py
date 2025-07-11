import math

def solve_orbital_time():
    """
    Calculates the time for a spacecraft to travel from the north pole
    to the equatorial plane of Mercury.
    """
    # Step 1: Define constants and convert units
    R_mercury_km = 2440  # km
    h_p_km = 218         # km
    h_a_km = 9982        # km
    mu_mercury_m3_s2 = 2.2032e13  # m^3/s^2

    R_m = R_mercury_km * 1000
    h_p_m = h_p_km * 1000
    h_a_m = h_a_km * 1000

    # Step 2: Calculate orbital parameters
    r_p = R_m + h_p_m
    r_a = R_m + h_a_m
    a = (r_p + r_a) / 2
    e = (r_a - r_p) / (r_a + r_p)

    # Step 3: Determine true anomalies for start and end points
    # From the problem description, we deduce:
    # Argument of periapsis, omega = 120 degrees
    # Start point (North Pole, lat=90deg) -> True anomaly, theta_1 = -30 degrees
    # End point (Equator, lat=0deg) -> True anomaly, theta_2 = 60 degrees
    theta1_rad = math.radians(-30)
    theta2_rad = math.radians(60)

    # Step 4: Calculate Time of Flight using Kepler's Equation
    
    # Function to convert true anomaly (theta) to eccentric anomaly (E)
    def true_to_eccentric(theta_rad, eccentricity):
        return math.atan2(math.sqrt(1 - eccentricity**2) * math.sin(theta_rad), eccentricity + math.cos(theta_rad))

    # Convert true anomalies to eccentric anomalies
    E1_rad = true_to_eccentric(theta1_rad, e)
    E2_rad = true_to_eccentric(theta2_rad, e)

    # Calculate mean anomalies
    M1_rad = E1_rad - e * math.sin(E1_rad)
    M2_rad = E2_rad - e * math.sin(E2_rad)
    
    # Calculate the mean motion (n)
    mean_motion_n = math.sqrt(mu_mercury_m3_s2 / a**3)
    
    # Calculate time of flight (delta_t)
    delta_t_s = (M2_rad - M1_rad) / mean_motion_n

    # Print the values used in the final calculation
    print("Values used in the final calculation:")
    print(f"Gravitational Parameter (μ): {mu_mercury_m3_s2} m^3/s^2")
    print(f"Semi-major axis (a): {a} m")
    print(f"Eccentricity (e): {e}")
    print(f"Start Eccentric Anomaly (E1): {E1_rad} rad")
    print(f"End Eccentric Anomaly (E2): {E2_rad} rad")

    # The final equation is: Time = ((E2 - e*sin(E2)) - (E1 - e*sin(E1))) / sqrt(μ / a^3)
    # We can write it as:
    # Time = (M2 - M1) / n
    # Time = (({E2_rad}) - ({e})*sin({E2_rad})) - (({E1_rad}) - ({e})*sin({E1_rad}))) / sqrt(({mu_mercury_m3_s2}) / ({a})**3)
    # Plugging in the numbers:
    print("\nFinal equation with values:")
    print(f"Time = (({E2_rad:.6f} - {e:.6f}*sin({E2_rad:.6f})) - ({E1_rad:.6f} - {e:.6f}*sin({E1_rad:.6f}))) / sqrt({mu_mercury_m3_s2} / {a}**3)")
    print(f"Time = {delta_t_s:.2f} seconds")

    # Step 5: Round the result to the nearest 10 seconds
    final_time = int(round(delta_t_s / 10.0) * 10)
    
    print("\nFinal Answer:")
    print(f"The calculated time is {delta_t_s:.2f} seconds.")
    print(f"Rounded to the nearest 10 seconds, the time is: {final_time}")
    
    return final_time

# Execute the function and print the final answer in the required format
final_answer = solve_orbital_time()
print(f'<<<{final_answer}>>>')
