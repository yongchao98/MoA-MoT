import math

def calculate_time_of_flight():
    """
    Calculates the time it took for a spacecraft to travel from the north pole
    through periapsis to the equatorial plane in orbit around Mercury.
    """
    # 1. Define constants and initial values in SI units.
    R_mercury = 2440 * 1000  # meters
    h_p = 218 * 1000         # meters, periapsis altitude
    h_a = 9982 * 1000        # meters, apoapsis altitude
    mu = 2.2032e13           # m^3 s^-2, standard gravitational parameter of Mercury

    # 2. Determine orbital elements.
    r_p = R_mercury + h_p  # Periapsis radius
    r_a = R_mercury + h_a  # Apoapsis radius
    a = (r_p + r_a) / 2    # Semi-major axis
    e = (r_a - r_p) / (r_a + r_p)  # Eccentricity

    # 3. Determine the true anomalies for the start and end points.
    # The path is 90N (North Pole) -> 60N (Periapsis) -> 0N (Equator).
    # This descending path means the argument of periapsis (omega) is 120 degrees.
    # We find the true anomaly (nu) for the start and end points using the relation
    # for a polar orbit: sin(latitude) = sin(omega + nu).
    # Start point (North Pole, 90N): nu_1 = -30 degrees
    # End point (Equator, 0N): nu_2 = 60 degrees
    nu1_deg = -30.0
    nu2_deg = 60.0

    # 4. Calculate time of flight using Kepler's Equation.

    # Convert true anomalies from degrees to radians for math functions.
    nu1_rad = math.radians(nu1_deg)
    nu2_rad = math.radians(nu2_deg)

    # Convert true anomalies (nu) to eccentric anomalies (E).
    # The relation is tan(E/2) = sqrt((1-e)/(1+e)) * tan(nu/2).
    sqrt_factor = math.sqrt((1 - e) / (1 + e))
    E1_rad = 2 * math.atan(sqrt_factor * math.tan(nu1_rad / 2))
    E2_rad = 2 * math.atan(sqrt_factor * math.tan(nu2_rad / 2))

    # Calculate time of flight (delta_t) using the eccentric anomalies.
    # Formula: Δt = (E₂ - E₁ - e * (sin(E₂) - sin(E₁))) * √(a³/μ)
    time_factor = math.sqrt((a**3) / mu)
    eccentric_anomaly_term = E2_rad - E1_rad - e * (math.sin(E2_rad) - math.sin(E1_rad))
    delta_t = eccentric_anomaly_term * time_factor

    # 5. Round the result to the nearest 10 seconds.
    final_time = int(round(delta_t / 10) * 10)

    # Print the explanation and the final equation with all numbers.
    print("The time of flight (Δt) is calculated using the formula:")
    print("Δt = (E₂ - E₁ - e * (sin(E₂) - sin(E₁))) * √(a³ / μ)\n")
    print("Here are the calculated values:")
    print(f"  Semi-major axis (a) = {a} m")
    print(f"  Eccentricity (e) = {e}")
    print(f"  Gravitational Parameter (μ) = {mu} m³/s²")
    print(f"  Initial Eccentric Anomaly (E₁) = {E1_rad} rad (from ν₁ = {nu1_deg}°)")
    print(f"  Final Eccentric Anomaly (E₂) = {E2_rad} rad (from ν₂ = {nu2_deg}°)")
    print("\nSubstituting the values into the equation:")
    print(f"Δt = ({E2_rad} - ({E1_rad}) - {e} * (sin({E2_rad}) - sin({E1_rad}))) * √({a}³ / {mu})")
    print(f"\nCalculated time of flight ≈ {delta_t:.2f} seconds")
    print(f"\nFinal Answer: Time rounded to the nearest 10 seconds is {final_time}")
    print(f"\n<<<{final_time}>>>")

# Run the calculation.
calculate_time_of_flight()