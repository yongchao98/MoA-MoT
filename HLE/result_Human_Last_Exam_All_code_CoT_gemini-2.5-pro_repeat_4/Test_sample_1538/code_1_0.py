import math

def solve_orbital_time():
    """
    Calculates the time for a spacecraft to travel from the North Pole
    through periapsis to the equatorial plane in its orbit around Mercury.
    """
    # Step 1: Establish the given physical constants and convert them to SI units.
    R_M = 2440 * 1000      # Radius of Mercury in meters
    mu = 2.2032e13         # Standard gravitational parameter of Mercury in m^3/s^2
    h_p = 218 * 1000       # Periapsis altitude in meters
    h_a = 9982 * 1000      # Apoapsis altitude in meters

    # Step 2: Calculate the orbital parameters.
    r_p = R_M + h_p
    r_a = R_M + h_a
    a = (r_p + r_a) / 2
    e = (r_a - r_p) / (r_a + r_p)
    n = math.sqrt(mu / a**3)

    # Step 3: Determine the true anomalies (ν) for the start and end points.
    # Start point: North Pole (latitude 90 deg), ν₁ = 30 deg.
    nu_start_rad = math.radians(30.0)
    # End point: Equatorial plane (latitude 0 deg), ν₂ = -60 deg.
    nu_end_rad = math.radians(-60.0)

    # Step 4: Calculate the time of flight using Kepler's laws.
    def calculate_time_from_periapsis(nu_rad, eccentricity, mean_motion):
        """
        Calculates time from periapsis passage for a given true anomaly.
        Returns the time (t), eccentric anomaly (E), and mean anomaly (M).
        """
        # Convert true anomaly (ν) to eccentric anomaly (E)
        tan_E_half = math.sqrt((1 - eccentricity) / (1 + eccentricity)) * math.tan(nu_rad / 2)
        E_rad = 2 * math.atan(tan_E_half)
        # Convert eccentric anomaly (E) to mean anomaly (M) using Kepler's Equation
        M_rad = E_rad - eccentricity * math.sin(E_rad)
        # Convert mean anomaly (M) to time (t)
        t = M_rad / mean_motion
        return t, E_rad, M_rad

    # Calculate time, eccentric anomaly, and mean anomaly for the start point
    t_start, E_start_rad, M_start_rad = calculate_time_from_periapsis(nu_start_rad, e, n)

    # Calculate time, eccentric anomaly, and mean anomaly for the end point
    t_end, E_end_rad, M_end_rad = calculate_time_from_periapsis(nu_end_rad, e, n)

    # The total travel time is the duration from North Pole to periapsis (t_start)
    # plus the duration from periapsis to the Equator (|t_end|).
    # This is equivalent to t_start - t_end, since t_end is negative.
    total_time = t_start - t_end

    # Step 5: Round the final answer to the nearest 10 seconds.
    final_answer = int(round(total_time / 10) * 10)

    # --- Output the results with each number in the final equations ---
    print("--- Orbital Parameter Calculation ---")
    print(f"Radius of periapsis (r_p) = {R_M:.0f} m + {h_p:.0f} m = {r_p:.0f} m")
    print(f"Radius of apoapsis (r_a) = {R_M:.0f} m + {h_a:.0f} m = {r_a:.0f} m")
    print(f"Semi-major axis (a) = ({r_p:.0f} m + {r_a:.0f} m) / 2 = {a:.0f} m")
    print(f"Eccentricity (e) = ({r_a:.0f} m - {r_p:.0f} m) / {r_p + r_a:.0f} m = {e:.6f}")
    print(f"Mean motion (n) = sqrt({mu:.4e} m^3/s^2 / {a**3:.4e} m^3) = {n:.6e} rad/s")

    print("\n--- Time Calculation for Start Point (North Pole) ---")
    print(f"True anomaly (ν₁) = {math.degrees(nu_start_rad):.1f} degrees")
    print(f"Eccentric anomaly (E₁) = {math.degrees(E_start_rad):.4f} degrees")
    print(f"Mean anomaly (M₁) = {E_start_rad:.6f} rad - {e:.6f} * sin({E_start_rad:.6f} rad) = {M_start_rad:.6f} rad")
    print(f"Time from periapsis (t₁) = {M_start_rad:.6f} rad / {n:.6e} rad/s = {t_start:.2f} s")

    print("\n--- Time Calculation for End Point (Equatorial Plane) ---")
    print(f"True anomaly (ν₂) = {math.degrees(nu_end_rad):.1f} degrees")
    print(f"Eccentric anomaly (E₂) = {math.degrees(E_end_rad):.4f} degrees")
    print(f"Mean anomaly (M₂) = {E_end_rad:.6f} rad - {e:.6f} * sin({E_end_rad:.6f} rad) = {M_end_rad:.6f} rad")
    print(f"Time from periapsis (t₂) = {M_end_rad:.6f} rad / {n:.6e} rad/s = {t_end:.2f} s")

    print("\n--- Total Time of Flight ---")
    print(f"Total Time = t₁ - t₂ = {t_start:.2f} s - ({t_end:.2f} s) = {total_time:.2f} s")

    print("\n--- Final Answer ---")
    print(f"Time rounded to the nearest 10 seconds is: {final_answer}")

solve_orbital_time()