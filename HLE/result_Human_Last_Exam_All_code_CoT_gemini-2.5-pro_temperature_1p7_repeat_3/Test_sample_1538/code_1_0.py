import math

def solve_flight_time():
    """
    Calculates the time of flight for a spacecraft in a polar orbit around Mercury.
    """
    # Step 0: Define constants
    R_mercury = 2440 * 1000  # Radius of Mercury in meters
    mu_mercury = 2.2032e13   # Standard gravitational parameter of Mercury in m^3 s^-2
    h_p_alt = 218 * 1000     # Periapsis altitude in meters
    h_a_alt = 9982 * 1000    # Apoapsis altitude in meters

    # Step 1: Calculate orbital parameters
    r_p = R_mercury + h_p_alt  # Radius of periapsis
    r_a = R_mercury + h_a_alt  # Radius of apoapsis
    a = (r_p + r_a) / 2      # Semi-major axis
    e = (r_a - r_p) / (r_a + r_p)  # Eccentricity

    # Step 2: Determine true anomalies for start and end points
    # Path: North Pole (90 deg N) -> Periapsis (60 deg N) -> Equator (0 deg N)
    # True anomaly (nu) is measured from periapsis (at 60 deg N).
    # Start point is 30 degrees before periapsis.
    nu1_deg = -30.0
    # End point is 60 degrees after periapsis.
    nu2_deg = 60.0
    
    nu1_rad = math.radians(nu1_deg)
    nu2_rad = math.radians(nu2_deg)

    # Step 3: Solve for time using Kepler's Equations
    
    # Function to get eccentric anomaly (E) from true anomaly (nu)
    def get_eccentric_anomaly(nu_rad, e):
        # Using the relation: cos(E) = (e + cos(nu)) / (1 + e * cos(nu))
        cos_E = (e + math.cos(nu_rad)) / (1 + e * math.cos(nu_rad))
        # acos gives E in [0, pi]. Adjust the sign to match the sign of nu.
        E_rad = math.acos(cos_E)
        if nu_rad < 0:
            E_rad = -E_rad
        return E_rad

    # Convert true anomalies to eccentric anomalies
    E1_rad = get_eccentric_anomaly(nu1_rad, e)
    E2_rad = get_eccentric_anomaly(nu2_rad, e)

    # Convert eccentric anomalies to mean anomalies using M = E - e*sin(E)
    M1_rad = E1_rad - e * math.sin(E1_rad)
    M2_rad = E2_rad - e * math.sin(E2_rad)

    # Calculate mean motion: n = sqrt(mu / a^3)
    n = math.sqrt(mu_mercury / a**3)

    # Calculate time of flight: delta_t = (M2 - M1) / n
    delta_t = (M2_rad - M1_rad) / n
    
    # Round the final answer to the nearest 10 seconds
    rounded_delta_t = int(round(delta_t / 10.0)) * 10

    # Output the results and the final calculation
    print("--- Orbital Parameters ---")
    print(f"Radius of periapsis (r_p): {r_p / 1000:.1f} km")
    print(f"Radius of apoapsis (r_a): {r_a / 1000:.1f} km")
    print(f"Semi-major axis (a): {a / 1000:.1f} km")
    print(f"Eccentricity (e): {e:.6f}")
    print("\n--- Time of Flight Calculation ---")
    print(f"Start True Anomaly (ν1): {nu1_deg}°")
    print(f"End True Anomaly (ν2): {nu2_deg}°")
    
    print("\n--- Final Equation: Time = (M2 - M1) / n ---")
    print(f"Mean Anomaly at start (M1): {M1_rad:.6f} rad")
    print(f"Mean Anomaly at end (M2): {M2_rad:.6f} rad")
    print(f"Mean Motion (n): {n:.9f} rad/s")
    print("\nPlugging in the numbers:")
    print(f"Time = ( {M2_rad:.6f} - ({M1_rad:.6f}) ) / {n:.9f}")
    print(f"Time = {M2_rad - M1_rad:.6f} / {n:.9f}")
    
    print(f"\nCalculated time of flight: {delta_t:.2f} seconds")
    print(f"Time rounded to the nearest 10 seconds: {rounded_delta_t} seconds")
    
    return rounded_delta_t

# Run the calculation and print the final answer in the required format
final_answer = solve_flight_time()
print(f"<<<{final_answer}>>>")
