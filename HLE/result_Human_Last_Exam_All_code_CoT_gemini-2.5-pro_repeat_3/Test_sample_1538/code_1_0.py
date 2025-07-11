import math

def solve_orbital_time():
    """
    Calculates the time for a spacecraft to travel from Mercury's north pole,
    through periapsis, to the equatorial plane.
    """
    # Step 1: Define constants and convert to SI units
    h_p_km = 218
    h_a_km = 9982
    R_M_km = 2440
    mu_m3_s2 = 2.2032e13

    h_p = h_p_km * 1000  # Periapsis altitude in meters
    h_a = h_a_km * 1000  # Apoapsis altitude in meters
    R_M = R_M_km * 1000  # Radius of Mercury in meters

    # Step 2: Calculate orbital parameters
    r_p = R_M + h_p
    r_a = R_M + h_a
    a = (r_p + r_a) / 2
    e = (r_a - r_p) / (r_a + r_p)
    n = math.sqrt(mu_m3_s2 / (a**3))

    # Step 3: Define true anomalies for start and end points
    # From the problem description, we deduced:
    # Start point (North Pole): nu_start = -30 degrees
    # End point (Equatorial Plane): nu_end = 60 degrees
    nu_start_deg = -30.0
    nu_end_deg = 60.0

    # Helper function to calculate time from periapsis for a given true anomaly
    def get_time_from_periapsis(nu_deg, eccentricity, mean_motion):
        """Calculates time from periapsis passage using Kepler's equation."""
        nu_rad = math.radians(nu_deg)
        
        # Calculate Eccentric Anomaly (E)
        E_half_tan = math.sqrt((1 - eccentricity) / (1 + eccentricity)) * math.tan(nu_rad / 2)
        E = 2 * math.atan(E_half_tan)
        
        # Calculate Mean Anomaly (M)
        M = E - eccentricity * math.sin(E)
        
        # Calculate time (t)
        t = M / mean_motion
        return t

    # Step 4: Calculate time of flight
    t_start = get_time_from_periapsis(nu_start_deg, e, n)
    t_end = get_time_from_periapsis(nu_end_deg, e, n)

    total_time = t_end - t_start

    # Step 5: Round and print the results
    rounded_time = int(round(total_time / 10) * 10)

    # Outputting the numbers in the final equation as requested
    print("--- Calculation Breakdown ---")
    print(f"Periapsis radius (r_p) = {R_M:.0f} m + {h_p:.0f} m = {r_p:.0f} m")
    print(f"Apoapsis radius (r_a) = {R_M:.0f} m + {h_a:.0f} m = {r_a:.0f} m")
    print(f"Semi-major axis (a) = ({r_p:.0f} m + {r_a:.0f} m) / 2 = {a:.0f} m")
    print(f"Eccentricity (e) = ({r_a:.0f} m - {r_p:.0f} m) / ({r_a:.0f} m + {r_p:.0f} m) = {e:.6f}")
    print(f"Mean motion (n) = sqrt({mu_m3_s2:.4e} m^3/s^2 / {a**3:.4e} m^3) = {n:.6e} rad/s\n")

    print(f"Time for spacecraft to reach North Pole (nu={nu_start_deg}°): {t_start:.1f} s from periapsis")
    print(f"Time for spacecraft to reach Equator (nu={nu_end_deg}°): {t_end:.1f} s from periapsis\n")

    print("--- Final Time Calculation ---")
    print(f"Total time = (Time to Equator) - (Time to North Pole)")
    print(f"Total time = {t_end:.1f} s - ({t_start:.1f} s) = {total_time:.1f} s\n")
    
    print(f"The final answer rounded to the nearest 10 seconds is: {rounded_time}")

solve_orbital_time()