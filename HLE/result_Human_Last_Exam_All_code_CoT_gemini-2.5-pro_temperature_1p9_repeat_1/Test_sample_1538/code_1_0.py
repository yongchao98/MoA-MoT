import math

def calculate_time_of_flight():
    """
    Calculates the time of flight for a spacecraft in a polar orbit around Mercury.
    The flight is from the North Pole (90 N), through periapsis (at 60 N), to the equatorial plane (0 N).
    """
    # Step 1: Define constants and given values
    R_M = 2440e3  # Mercury's radius in meters
    mu_M = 2.2032e13  # Mercury's standard gravitational parameter in m^3/s^2
    h_p = 218e3  # Periapsis altitude in meters
    h_a = 9982e3  # Apoapsis altitude in meters

    # Step 2: Calculate orbital parameters
    r_p = R_M + h_p  # Periapsis radius
    r_a = R_M + h_a  # Apoapsis radius
    a = (r_p + r_a) / 2  # Semi-major axis
    e = (r_a - r_p) / (r_a + r_p)  # Eccentricity
    n = math.sqrt(mu_M / a**3)  # Mean motion

    # Step 3: Determine true anomalies (nu) in radians
    # Initial point: North Pole (90 N), which is 30 degrees before periapsis (60 N)
    nu1_rad = math.radians(-30)
    # Final point: Equatorial Plane (0 N), which is 60 degrees after periapsis (60 N)
    nu2_rad = math.radians(60)

    # Step 4: Solve for Time of Flight using Kepler's Equation
    
    # Function to calculate eccentric anomaly (E) from true anomaly (nu)
    def get_eccentric_anomaly(nu, e):
        tan_nu_half = math.tan(nu / 2)
        tan_E_half = math.sqrt((1 - e) / (1 + e)) * tan_nu_half
        return 2 * math.atan(tan_E_half)

    # Calculate eccentric anomalies for start and end points
    E1 = get_eccentric_anomaly(nu1_rad, e)
    E2 = get_eccentric_anomaly(nu2_rad, e)

    # Calculate mean anomalies for start and end points using M = E - e*sin(E)
    M1 = E1 - e * math.sin(E1)
    M2 = E2 - e * math.sin(E2)

    # The final equation for time of flight is delta_t = (M2 - M1) / n
    # We will now print each value in this equation
    print("Values for the time of flight equation: delta_t = (M2 - M1) / n")
    print(f"Final Mean Anomaly (M2): {M2} radians")
    print(f"Initial Mean Anomaly (M1): {M1} radians")
    print(f"Mean Motion (n): {n} radians/s")
    
    # Calculate total time of flight
    delta_t = (M2 - M1) / n

    # Step 5: Format and print the final answer
    # Round to the nearest 10 seconds
    final_answer = int(round(delta_t / 10) * 10)

    print("\n--- Final Calculation ---")
    print(f"The total time of flight is ({M2}) - ({M1})) / {n} = {delta_t:.2f} seconds.")
    print(f"Rounding to the nearest 10 seconds, the time taken is: {final_answer}")
    
    return final_answer

# Execute the calculation and print the final answer in the required format
final_time = calculate_time_of_flight()
print(f"<<<{final_time}>>>")
