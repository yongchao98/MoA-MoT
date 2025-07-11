import math

def solve_spacecraft_time():
    """
    Calculates the time it took for a spacecraft to travel from the north pole
    through periapsis to the equatorial plane around Mercury.
    """
    # --- Given constants ---
    # Radius of Mercury in meters
    R_M = 2440 * 1000
    # Standard gravitational parameter of Mercury in m^3/s^2
    mu = 2.2032e13
    # Periapsis altitude in meters
    h_p = 218 * 1000
    # Apoapsis altitude in meters
    h_a = 9982 * 1000

    # --- Step 1: Calculate Orbital Parameters ---
    # Radii of periapsis and apoapsis from the center of Mercury
    r_p = R_M + h_p
    r_a = R_M + h_a

    # Semi-major axis (a)
    a = (r_p + r_a) / 2
    # Eccentricity (e)
    e = (r_a - r_p) / (r_a + r_p)

    # --- Step 2 & 3: True Anomalies ---
    # As determined by geometric analysis of the path "from North Pole,
    # through periapsis, to equator" for a polar orbit with periapsis at 60N.
    # Start true anomaly (nu1) at the North Pole
    nu1_deg = -30
    # End true anomaly (nu2) at the Equator
    nu2_deg = 60
    
    nu1_rad = math.radians(nu1_deg)
    nu2_rad = math.radians(nu2_deg)

    # --- Step 4: Calculate Time of Flight using Kepler's Equation ---

    # Function to convert true anomaly (nu) to eccentric anomaly (E)
    def true_to_eccentric_anomaly(nu, e):
        # The relationship is tan(E/2) = sqrt((1-e)/(1+e)) * tan(nu/2)
        tan_E_half = math.sqrt((1 - e) / (1 + e)) * math.tan(nu / 2)
        return 2 * math.atan(tan_E_half)

    # Function to convert eccentric anomaly (E) to mean anomaly (M)
    def eccentric_to_mean_anomaly(E, e):
        # Kepler's Equation: M = E - e*sin(E)
        return E - e * math.sin(E)

    # Calculate eccentric anomalies E1 and E2
    E1 = true_to_eccentric_anomaly(nu1_rad, e)
    E2 = true_to_eccentric_anomaly(nu2_rad, e)

    # Calculate mean anomalies M1 and M2
    M1 = eccentric_to_mean_anomaly(E1, e)
    M2 = eccentric_to_mean_anomaly(E2, e)
    
    # Calculate the change in mean anomaly
    delta_M = M2 - M1
    
    # Calculate the mean motion (n)
    n = math.sqrt(mu / a**3)
    
    # Calculate the time of flight (delta_t)
    delta_t = delta_M / n
    
    # Round the final answer to the nearest 10 seconds
    final_time = int(round(delta_t, -1))

    # --- Output the results ---
    print("Calculation Details:")
    print(f"Semi-major axis (a): {a} m")
    print(f"Eccentricity (e): {e}")
    print(f"Mean motion (n): {n} rad/s")
    print(f"Change in mean anomaly (ΔM): {delta_M} rad")
    print(f"\nTime of flight calculation:")
    print(f"Δt = ΔM / n")
    print(f"Δt = {delta_M} / {n}")
    print(f"Δt = {delta_t:.2f} seconds")
    print(f"\nFinal answer rounded to the nearest 10 seconds: {final_time}")
    
    return final_time

# Run the calculation and store the final answer
final_answer = solve_spacecraft_time()
print(f"<<<{final_answer}>>>")