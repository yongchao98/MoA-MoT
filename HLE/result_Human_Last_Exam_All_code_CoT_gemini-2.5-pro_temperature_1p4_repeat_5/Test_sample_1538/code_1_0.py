import math

def solve_time_of_flight():
    """
    Calculates the time it took for the spacecraft to travel from the north pole
    through periapsis to the equatorial plane.
    """
    # Step 1: Define constants and convert to SI units
    h_p_km = 218         # Periapsis altitude in km
    h_a_km = 9982        # Apoapsis altitude in km
    R_km = 2440          # Radius of Mercury in km
    mu_m3_s2 = 2.2032e13 # Gravitational parameter of Mercury in m^3/s^-2

    h_p_m = h_p_km * 1000
    h_a_m = h_a_km * 1000
    R_m = R_km * 1000

    print("### Solving for Spacecraft Time of Flight ###\n")
    print("--- Given Parameters ---")
    print(f"Periapsis Altitude (h_p): {h_p_km} km")
    print(f"Apoapsis Altitude (h_a): {h_a_km} km")
    print(f"Mercury Radius (R): {R_km} km")
    print(f"Mercury Gravitational Parameter (μ): {mu_m3_s2} m^3 s^-2\n")

    # Step 2: Calculate orbital parameters
    r_p = R_m + h_p_m
    r_a = R_m + h_a_m
    a = (r_p + r_a) / 2
    e = (r_a - r_p) / (r_a + r_p)
    n = math.sqrt(mu_m3_s2 / (a**3))
    
    print("--- Calculated Orbital Parameters ---")
    print(f"Semi-major axis (a): {a:.3f} m")
    print(f"Eccentricity (e): {e:.6f}")
    print(f"Mean motion (n): {n:.9f} rad/s\n")

    # Step 3: Define true anomalies for start and end points
    nu1_deg = 30.0   # True anomaly at North Pole
    nu2_deg = -60.0  # True anomaly at Equator crossing

    def get_time_from_periapsis(nu_deg, e, n):
        """Calculates time from periapsis for a given true anomaly."""
        nu_rad = math.radians(nu_deg)
        
        # Eccentric anomaly E from true anomaly nu
        E_rad = 2 * math.atan2(math.sqrt(1 - e) * math.sin(nu_rad / 2),
                               math.sqrt(1 + e) * math.cos(nu_rad / 2))
        
        # Mean anomaly M from eccentric anomaly E (Kepler's Equation)
        M_rad = E_rad - e * math.sin(E_rad)
        
        # Time t from mean anomaly M
        time = M_rad / n
        return time

    # Step 4: Calculate time for each segment
    
    # Time from periapsis to North Pole (nu=30 deg)
    time_segment1 = get_time_from_periapsis(nu1_deg, e, n)
    
    # Time from periapsis to Equator (nu=-60 deg)
    time_segment2 = get_time_from_periapsis(nu2_deg, e, n)

    # Durations are always positive
    duration1 = abs(time_segment1)
    duration2 = abs(time_segment2)
    
    # Total time is the sum of the durations of the two segments
    total_time_sec = duration1 + duration2
    
    print("--- Time of Flight Calculation ---")
    print("The total time is the sum of two segments:")
    print("1. Time from North Pole (ν = 30°) to Periapsis (ν = 0°)")
    print("2. Time from Periapsis (ν = 0°) to Equator (ν = -60°)\n")

    print(f"Time for segment 1: {duration1:.2f} s")
    print(f"Time for segment 2: {duration2:.2f} s")
    print("\n--- Final Equation ---")
    print(f"Total Time = {duration1:.2f} + {duration2:.2f}")
    
    # Step 5: Final calculation and rounding
    final_answer = int(round(total_time_sec, -1))
    
    print(f"\nResult: {total_time_sec:.2f} seconds")
    print(f"Rounded to the nearest 10 seconds, the final answer is: {final_answer}")

solve_time_of_flight()
<<<1270>>>