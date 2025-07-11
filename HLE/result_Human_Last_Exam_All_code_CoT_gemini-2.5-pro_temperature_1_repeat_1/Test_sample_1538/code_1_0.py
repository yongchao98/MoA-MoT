import math

def solve_orbital_time():
    """
    Calculates the time of flight for a spacecraft in a polar orbit around Mercury
    from the North Pole to the Equator, passing through periapsis.
    """
    # --- Given Constants ---
    # Altitudes and Radius in meters
    h_p = 218 * 1000
    h_a = 9982 * 1000
    R_M = 2440 * 1000
    # Standard gravitational parameter of Mercury in m^3/s^2
    mu = 2.2032e13

    # --- Step 1: Calculate Orbital Parameters ---
    # Radius of periapsis and apoapsis
    r_p = R_M + h_p
    r_a = R_M + h_a
    # Semi-major axis (a)
    a = (r_p + r_a) / 2
    # Eccentricity (e)
    e = (r_a - r_p) / (r_a + r_p)
    # Mean motion (n) in rad/s
    n = math.sqrt(mu / a**3)

    # --- Step 2: Determine Orbital Orientation and True Anomalies ---
    # For a polar orbit (i=90 deg), sin(latitude) = sin(argument_of_periapsis + true_anomaly).
    # Since the spacecraft is moving from a higher latitude (90N) to a lower one (60N),
    # the argument of periapsis (omega) must be 120 degrees for the latitude to be decreasing at periapsis.
    omega_deg = 120
    
    # At the start point (North Pole, latitude = 90 deg):
    # sin(90) = sin(120 + nu_start) => 120 + nu_start = 90 => nu_start = -30 deg
    nu_start_deg = -30
    
    # At the end point (Equator, latitude = 0 deg):
    # sin(0) = sin(120 + nu_end) => 120 + nu_end = 180 => nu_end = 60 deg
    nu_end_deg = 60
    
    # Convert true anomalies to radians for calculations
    nu_start_rad = math.radians(nu_start_deg)
    nu_end_rad = math.radians(nu_end_deg)

    # --- Step 3: Convert True Anomalies (nu) to Eccentric (E) and Mean (M) Anomalies ---
    
    # Helper function to find Eccentric Anomaly E from True Anomaly nu
    def get_eccentric_anomaly(nu_rad, e):
        tan_nu_half = math.tan(nu_rad / 2)
        tan_E_half = math.sqrt((1 - e) / (1 + e)) * tan_nu_half
        return 2 * math.atan(tan_E_half)

    # Helper function to find Mean Anomaly M from Eccentric Anomaly E using Kepler's Equation
    def get_mean_anomaly(E_rad, e):
        return E_rad - e * math.sin(E_rad)

    # Calculate E and M for the start and end points
    E_start_rad = get_eccentric_anomaly(nu_start_rad, e)
    M_start_rad = get_mean_anomaly(E_start_rad, e)

    E_end_rad = get_eccentric_anomaly(nu_end_rad, e)
    M_end_rad = get_mean_anomaly(E_end_rad, e)

    # --- Step 4: Calculate Time of Flight (Δt) ---
    # The final equation is Δt = (M_end - M_start) / n
    delta_t_seconds = (M_end_rad - M_start_rad) / n

    # --- Step 5: Round and Print the Final Answer ---
    # Round the result to the nearest 10 seconds
    final_time = int(round(delta_t_seconds / 10) * 10)

    print("Calculation of Time of Flight (Δt)")
    print("-" * 40)
    print("The final equation is: Δt = (M_end - M_start) / n")
    print("\nValues used in the final equation:")
    print(f"  Mean Anomaly at End Point (M_end): {M_end_rad:.8f} rad")
    print(f"  Mean Anomaly at Start Point (M_start): {M_start_rad:.8f} rad")
    print(f"  Mean Motion (n): {n:.8f} rad/s")

    print("\nSubstituting the numbers into the equation:")
    print(f"  Δt = ({M_end_rad:.8f} - ({M_start_rad:.8f})) / {n:.8f}")
    print(f"  Δt = {(M_end_rad - M_start_rad):.8f} / {n:.8f}")
    print(f"  Δt = {delta_t_seconds:.2f} seconds")

    print("\nFinal Answer")
    print("-" * 40)
    print(f"The calculated time of flight, rounded to the nearest 10 seconds, is: {final_time}")
    
    return final_time

# Execute the function and store the result for the final formatted output
final_answer = solve_orbital_time()
print(f"<<<{final_answer}>>>")
