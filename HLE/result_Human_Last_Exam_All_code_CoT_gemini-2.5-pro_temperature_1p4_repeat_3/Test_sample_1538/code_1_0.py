import numpy as np

def solve_orbital_time():
    """
    Calculates the time for a spacecraft to travel from the North Pole,
    through periapsis, to the equatorial plane of Mercury.
    """
    # Step 1: Define constants and convert units
    h_p_km = 218
    h_a_km = 9982
    R_km = 2440
    mu_m3_s2 = 2.2032e13

    # Convert all length units to meters
    h_p = h_p_km * 1000
    h_a = h_a_km * 1000
    R = R_km * 1000

    # Step 2: Calculate orbital parameters
    r_p = R + h_p  # Radius of periapsis
    r_a = R + h_a  # Radius of apoapsis
    a = (r_p + r_a) / 2  # Semi-major axis
    e = (r_a - r_p) / (r_a + r_p)  # Eccentricity

    # Step 3: Determine true anomalies for the start and end points
    # The path is North Pole (90N) -> Periapsis (60N) -> Equator (0).
    # From the formula sin(latitude) = sin(argument_of_periapsis + true_anomaly)
    # and the path description, we deduce the start and end true anomalies.
    nu1_deg = -30  # True anomaly at the North Pole
    nu2_deg = 60   # True anomaly at the Equator (descending node)

    # Convert degrees to radians for calculations
    nu1_rad = np.radians(nu1_deg)
    nu2_rad = np.radians(nu2_deg)

    # Step 4: Calculate time of flight using Kepler's Laws
    
    # Define a helper function to get mean anomaly from true anomaly
    def get_mean_anomaly(nu_rad, e):
        # Calculate eccentric anomaly (E) from true anomaly (nu)
        E = 2 * np.arctan(np.sqrt((1 - e) / (1 + e)) * np.tan(nu_rad / 2))
        # Calculate mean anomaly (M) from eccentric anomaly (E)
        M = E - e * np.sin(E)
        return E, M

    # Calculate eccentric and mean anomalies for the start and end points
    E1_rad, M1_rad = get_mean_anomaly(nu1_rad, e)
    E2_rad, M2_rad = get_mean_anomaly(nu2_rad, e)

    # The time of flight (delta_t) is the change in mean anomaly divided by mean motion (n)
    # delta_t = (M2 - M1) / n  where n = sqrt(mu / a^3)
    delta_M = M2_rad - M1_rad
    time_of_flight = delta_M * np.sqrt(a**3 / mu_m3_s2)

    # Step 5: Finalize and print the result
    # Round the result to the nearest 10 seconds
    final_time = int(round(time_of_flight / 10.0)) * 10
    
    # Print the detailed calculation as requested
    print(f"Calculation for time of flight from North Pole to Equator via Periapsis")
    print("-" * 60)
    print(f"Given Parameters:")
    print(f"  Periapsis altitude (h_p) = {h_p_km} km")
    print(f"  Apoapsis altitude (h_a) = {h_a_km} km")
    print(f"  Mercury radius (R) = {R_km} km")
    print(f"  Mercury gravitational parameter (mu) = {mu_m3_s2} m^3/s^2\n")
    
    print(f"Derived Orbital Parameters:")
    print(f"  Semi-major axis (a) = {a} m")
    print(f"  Eccentricity (e) = {e:.6f}\n")

    print(f"Anomalies and Time Calculation:")
    print(f"  Start point (North Pole): True Anomaly (nu1) = {nu1_deg} degrees")
    print(f"    -> Eccentric Anomaly (E1) = {E1_rad:.6f} rad")
    print(f"    -> Mean Anomaly (M1) = {M1_rad:.6f} rad\n")
    
    print(f"  End point (Equator): True Anomaly (nu2) = {nu2_deg} degrees")
    print(f"    -> Eccentric Anomaly (E2) = {E2_rad:.6f} rad")
    print(f"    -> Mean Anomaly (M2) = {M2_rad:.6f} rad\n")

    print(f"Final Equation:")
    print(f"Time = (M2 - M1) * sqrt(a^3 / mu)")
    print(f"Time = ({M2_rad:.6f} - ({M1_rad:.6f})) * sqrt({a}^3 / {mu_m3_s2})")
    print(f"Time = {time_of_flight:.2f} seconds\n")
    
    print(f"Result rounded to the nearest 10 seconds: {final_time}")
    
    return final_time

# Execute the function and print the final answer in the required format
final_answer = solve_orbital_time()
print(f"<<<{final_answer}>>>")
