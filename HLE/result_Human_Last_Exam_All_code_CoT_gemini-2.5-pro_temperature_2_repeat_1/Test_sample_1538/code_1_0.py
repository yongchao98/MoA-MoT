import math

def solve_orbital_time():
    """
    Calculates the time for a spacecraft to travel from the North Pole,
    through periapsis, to the equatorial plane of Mercury.
    """
    # Step 1: Define given constants and convert to SI units (meters)
    R_mercury_km = 2440
    mu_mercury_m3_s2 = 2.2032e13
    h_p_km = 218
    h_a_km = 9982

    R = R_mercury_km * 1000  # Radius of Mercury in meters
    h_p = h_p_km * 1000     # Periapsis altitude in meters
    h_a = h_a_km * 1000     # Apoapsis altitude in meters

    # Step 2: Calculate orbital parameters
    r_p = R + h_p
    r_a = R + h_a
    a = (r_p + r_a) / 2  # Semi-major axis
    e = (r_a - r_p) / (r_a + r_p) # Eccentricity

    # Step 3: Define start and end true anomalies in radians
    # As explained in the plan, the journey is from the North Pole (theta_1 = -30 deg)
    # to the equatorial plane (theta_2 = 60 deg) after passing periapsis (theta = 0 deg).
    theta1_deg = -30
    theta2_deg = 60
    theta1 = math.radians(theta1_deg)
    theta2 = math.radians(theta2_deg)
    
    # Step 4: Calculate time of flight

    # Function to calculate eccentric anomaly (E) from true anomaly (theta)
    def get_eccentric_anomaly(theta, e):
        # Using atan2 to ensure the correct quadrant for the angle
        cos_E = (e + math.cos(theta)) / (1 + e * math.cos(theta))
        sin_E = (math.sqrt(1 - e**2) * math.sin(theta)) / (1 + e * math.cos(theta))
        return math.atan2(sin_E, cos_E)

    # Calculate eccentric anomalies for start and end points
    E1 = get_eccentric_anomaly(theta1, e)
    E2 = get_eccentric_anomaly(theta2, e)

    # Calculate mean anomalies using Kepler's Equation: M = E - e*sin(E)
    M1 = E1 - e * math.sin(E1)
    M2 = E2 - e * math.sin(E2)
    
    # Calculate mean motion (n)
    n = math.sqrt(mu_mercury_m3_s2 / a**3)

    # Calculate the time of flight (delta_t)
    delta_t = (M2 - M1) / n

    # Step 5: Round the result to the nearest 10 seconds
    final_time = int(round(delta_t / 10.0) * 10)

    print("--- Orbital Parameters ---")
    print(f"Radius of Mercury: {R_mercury_km} km")
    print(f"Periapsis altitude (h_p): {h_p_km} km")
    print(f"Apoapsis altitude (h_a): {h_a_km} km")
    print(f"Standard gravitational parameter (μ): {mu_mercury_m3_s2:.4e} m^3 s^-2")
    print(f"Radius of periapsis (r_p): {r_p:.0f} m")
    print(f"Radius of apoapsis (r_a): {r_a:.0f} m")
    print(f"Semi-major axis (a): {a:.0f} m")
    print(f"Eccentricity (e): {e:.6f}")
    print("\n--- Time of Flight Calculation ---")
    print(f"Start true anomaly (θ₁): {theta1_deg}°")
    print(f"End true anomaly (θ₂): {theta2_deg}°")
    print(f"Start eccentric anomaly (E₁): {math.degrees(E1):.4f}°")
    print(f"End eccentric anomaly (E₂): {math.degrees(E2):.4f}°")
    print(f"Start mean anomaly (M₁): {M1:.6f} rad")
    print(f"End mean anomaly (M₂): {M2:.6f} rad")
    print(f"Mean motion (n): {n:.6e} rad/s")
    print(f"Time of flight (unrounded): {delta_t:.4f} s")
    print("\n--- Final Answer ---")
    print(f"The time it took for the spacecraft to travel is {final_time} seconds (rounded to the nearest 10).")
    
    # Return the final integer answer for the "answer" block
    return final_time

# Run the calculation and store the result
calculated_answer = solve_orbital_time()

# Final answer in the required format
# print(f"<<<{calculated_answer}>>>")