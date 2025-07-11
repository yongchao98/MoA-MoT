import math

def solve_orbit_time():
    """
    Calculates the time for a spacecraft to travel from the North Pole,
    through periapsis, to the equatorial plane of Mercury.
    """
    # --- Step 0: Define Given Constants ---
    # Altitudes converted to meters
    h_p = 218 * 1000
    # Apoapsis altitude converted to meters
    h_a = 9982 * 1000
    # Radius of Mercury in meters
    R_merc = 2440 * 1000
    # Standard gravitational parameter of Mercury in m^3/s^2
    mu_merc = 2.2032e13

    print("--- Input Parameters ---")
    print(f"Periapsis altitude (h_p): {h_p / 1000} km")
    print(f"Apoapsis altitude (h_a): {h_a / 1000} km")
    print(f"Radius of Mercury (R): {R_merc / 1000} km")
    print(f"Gravitational Parameter (μ): {mu_merc} m^3/s^2\n")

    # --- Step 1: Calculate Orbital Parameters ---
    # Radii of periapsis and apoapsis
    r_p = R_merc + h_p
    r_a = R_merc + h_a

    # Semi-major axis (a)
    a = (r_p + r_a) / 2

    # Eccentricity (e)
    e = (r_a - r_p) / (r_a + r_p)

    # Mean motion (n) in rad/s
    n = math.sqrt(mu_merc / a**3)
    
    print("--- Calculated Orbital Parameters ---")
    print(f"Radius of periapsis (r_p) = {R_merc} + {h_p} = {r_p} m")
    print(f"Radius of apoapsis (r_a) = {R_merc} + {h_a} = {r_a} m")
    print(f"Semi-major axis (a) = ({r_p} + {r_a}) / 2 = {a} m")
    print(f"Eccentricity (e) = ({r_a} - {r_p}) / ({r_p} + {r_a}) = {e}")
    print(f"Mean motion (n) = sqrt({mu_merc} / {a}^3) = {n} rad/s\n")

    # --- Step 2: Determine True Anomalies (ν) ---
    # As derived in the plan, the journey starts at the North Pole (ν1 = -30°)
    # and ends at the Equator (ν2 = 60°).
    nu1_deg = -30.0
    nu2_deg = 60.0
    
    # Convert true anomalies to radians for calculation
    nu1_rad = math.radians(nu1_deg)
    nu2_rad = math.radians(nu2_deg)
    
    print("--- True Anomalies (ν) ---")
    print(f"Start true anomaly (ν1): {nu1_deg} degrees")
    print(f"End true anomaly (ν2): {nu2_deg} degrees\n")

    # --- Step 3: Calculate Eccentric Anomalies (E) ---
    # Relationship: tan(E/2) = sqrt((1-e)/(1+e)) * tan(ν/2)
    conv_factor = math.sqrt((1 - e) / (1 + e))
    
    # For starting point (ν1)
    tan_E1_over_2 = conv_factor * math.tan(nu1_rad / 2)
    E1_rad = 2 * math.atan(tan_E1_over_2)

    # For ending point (ν2)
    tan_E2_over_2 = conv_factor * math.tan(nu2_rad / 2)
    E2_rad = 2 * math.atan(tan_E2_over_2)
    
    print("--- Eccentric Anomalies (E) ---")
    print(f"Start eccentric anomaly (E1): {math.degrees(E1_rad)} degrees")
    print(f"End eccentric anomaly (E2): {math.degrees(E2_rad)} degrees\n")

    # --- Step 4: Calculate Mean Anomalies (M) ---
    # Kepler's Equation: M = E - e*sin(E)
    M1_rad = E1_rad - e * math.sin(E1_rad)
    M2_rad = E2_rad - e * math.sin(E2_rad)
    
    print("--- Mean Anomalies (M) ---")
    print(f"Start mean anomaly (M1) = {E1_rad} - {e}*sin({E1_rad}) = {M1_rad} rad")
    print(f"End mean anomaly (M2) = {E2_rad} - {e}*sin({E2_rad}) = {M2_rad} rad\n")

    # --- Step 5: Calculate Time of Flight ---
    delta_M = M2_rad - M1_rad
    time_of_flight_s = delta_M / n
    
    print("--- Time of Flight Calculation ---")
    print(f"Time = (M2 - M1) / n")
    print(f"Time = ({M2_rad} - ({M1_rad})) / {n}")
    print(f"Time = {delta_M} / {n} = {time_of_flight_s} s\n")

    # --- Step 6: Round and Output Final Answer ---
    rounded_time = int(round(time_of_flight_s / 10) * 10)
    
    print("--- Final Answer ---")
    print("The time it took for the spacecraft to travel from the north pole through periapsis to the equatorial plane is:")
    print(rounded_time)

solve_orbit_time()
<<<1270>>>