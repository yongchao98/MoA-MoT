import math

def solve_orbit_time():
    """
    Calculates the time of flight for a spacecraft in orbit around Mercury.
    """
    # Step 1: Define constants and convert units to meters and seconds.
    R_MERCURY = 2440 * 1000  # Radius of Mercury in meters
    MU_MERCURY = 2.2032e13   # Standard gravitational parameter of Mercury in m^3/s^2
    h_p = 218 * 1000         # Periapsis altitude in meters
    h_a = 9982 * 1000        # Apoapsis altitude in meters

    # Step 2: Calculate orbital parameters.
    # Radius of periapsis and apoapsis
    r_p = R_MERCURY + h_p
    r_a = R_MERCURY + h_a

    # Semi-major axis
    a = (r_p + r_a) / 2

    # Eccentricity
    e = (r_a - r_p) / (r_a + r_p)

    # Mean motion
    n = math.sqrt(MU_MERCURY / a**3)

    # Step 3: Define start and end true anomalies based on the geometric path.
    # As explained in the plan, the path from North Pole -> Periapsis -> Equator corresponds
    # to a flight from a true anomaly of -30 degrees to +60 degrees.
    nu_start_deg = -30.0
    nu_end_deg = 60.0

    nu_start_rad = math.radians(nu_start_deg)
    nu_end_rad = math.radians(nu_end_deg)

    # Helper function to calculate time from periapsis for a given true anomaly
    def get_time_from_periapsis(nu_rad, a_val, e_val, n_val):
        # Calculate eccentric anomaly (E)
        tan_E_half = math.sqrt((1 - e_val) / (1 + e_val)) * math.tan(nu_rad / 2)
        E_rad = 2 * math.atan(tan_E_half)
        # Calculate mean anomaly (M)
        M_rad = E_rad - e_val * math.sin(E_rad)
        # Calculate time (t)
        t = M_rad / n_val
        return t, E_rad, M_rad

    # Step 4: Calculate the time from periapsis for the start and end points.
    t_start, E_start, M_start = get_time_from_periapsis(nu_start_rad, a, e, n)
    t_end, E_end, M_end = get_time_from_periapsis(nu_end_rad, a, e, n)

    # Step 5: Calculate the total time of flight.
    total_time = t_end - t_start

    # Step 6: Round the final answer to the nearest 10 seconds.
    rounded_time = int(round(total_time / 10) * 10)

    # Step 7: Output the calculation steps and the final answer.
    print("The time of flight is calculated as Δt = t_end - t_start, where t is the time from periapsis.\n")
    print(f"Key orbital parameters: a = {a:.0f} m, e = {e:.6f}, n = {n:.8f} rad/s\n")
    
    print(f"Time calculation for start point (North Pole, ν = {nu_start_deg}°):")
    print(f"t_start = (E_start - e*sin(E_start)) / n")
    print(f"t_start = ({E_start:.6f} - {e:.6f}*sin({E_start:.6f})) / {n:.8f}")
    print(f"t_start = {M_start:.6f} / {n:.8f} = {t_start:.2f} s\n")

    print(f"Time calculation for end point (Equator, ν = {nu_end_deg}°):")
    print(f"t_end = (E_end - e*sin(E_end)) / n")
    print(f"t_end = ({E_end:.6f} - {e:.6f}*sin({E_end:.6f})) / {n:.8f}")
    print(f"t_end = {M_end:.6f} / {n:.8f} = {t_end:.2f} s\n")
    
    print("Total time of flight (Δt):")
    print(f"Δt = {t_end:.2f} - ({t_start:.2f})")
    print(f"Δt = {total_time:.2f} s\n")

    print(f"The final answer rounded to the nearest 10 seconds is: {rounded_time}")

    # Final answer in the required format
    print(f"\n<<<{rounded_time}>>>")

solve_orbit_time()