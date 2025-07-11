import math

def solve_orbit_time():
    """
    Calculates the time for a spacecraft to travel from the North Pole 
    to the equator through periapsis in a specific orbit around Mercury.
    """
    # --- GIVENS ---
    # Convert all initial values to SI units (meters, seconds)
    h_p = 218 * 1000  # Periapsis altitude in meters
    h_a = 9982 * 1000 # Apoapsis altitude in meters
    R = 2440 * 1000   # Radius of Mercury in meters
    mu = 2.2032e13    # Standard gravitational parameter of Mercury in m^3/s^2

    # As determined by the problem geometry:
    nu1_deg = 30.0   # True anomaly at North Pole (start point)
    nu2_deg = -60.0  # True anomaly at Equator (end point)

    # --- 1. CALCULATE ORBITAL PARAMETERS ---
    r_p = R + h_p  # Periapsis radius
    r_a = R + h_a  # Apoapsis radius
    a = (r_p + r_a) / 2.0  # Semi-major axis
    e = (r_a - r_p) / (r_a + r_p) # Eccentricity
    n = math.sqrt(mu / a**3)  # Mean motion in rad/s

    # --- 2. CALCULATE MEAN ANOMALIES (M1, M2) ---
    def get_mean_anomaly(nu_deg, eccentricity):
        """Calculates mean anomaly from true anomaly and eccentricity."""
        nu_rad = math.radians(nu_deg)
        # Relation between Eccentric (E) and True (nu) anomaly
        tan_E_div_2 = math.sqrt((1 - eccentricity) / (1 + eccentricity)) * math.tan(nu_rad / 2)
        E_rad = 2 * math.atan(tan_E_div_2)
        # Kepler's Equation: Mean Anomaly M = E - e*sin(E)
        M_rad = E_rad - eccentricity * math.sin(E_rad)
        return M_rad

    M1_rad = get_mean_anomaly(nu1_deg, e)
    M2_rad = get_mean_anomaly(nu2_deg, e)
    
    # --- 3. CALCULATE TIME OF FLIGHT ---
    # The time of flight is the change in mean anomaly divided by the mean motion.
    time_seconds = (M1_rad - M2_rad) / n

    # --- 4. ROUND and PREPARE OUTPUT ---
    # Round to the nearest 10 seconds
    rounded_time = int(round(time_seconds / 10.0) * 10.0)
    
    # --- 5. PRINT THE RESULTS AND FINAL EQUATION ---
    print("This script calculates the time of flight for a spacecraft orbiting Mercury.")
    print("All calculations use SI units and no intermediate rounding is performed.\n")
    print("--- Calculated Orbital Parameters ---")
    print(f"Semi-major axis (a): {a} m")
    print(f"Eccentricity (e): {e}")
    print(f"Mean motion (n): {n} rad/s\n")

    print("--- Time of Flight Calculation ---")
    print(f"The calculation for the time of flight (t) is based on the equation: t = (M1 - M2) / n")
    print(f"Where M1 is the mean anomaly at the start point (nu={nu1_deg} deg) and M2 is the mean anomaly at the end point (nu={nu2_deg} deg).")
    
    print("\nThe final equation with the calculated values is:")
    print(f"t = ({M1_rad} - ({M2_rad})) / {n}")
    
    print(f"\nThe unrounded time of flight is {time_seconds} seconds.")
    print(f"The final answer, rounded to the nearest 10 seconds, is {rounded_time} seconds.")
    print(f"\n<<<{rounded_time}>>>")


solve_orbit_time()