import math

def calculate_incident_power():
    """
    Calculates the total power incident on a photovoltaic cell on the Moon
    after light from a source is reflected by two satellites.
    """
    # Given data
    P = 1e9  # Source Luminosity in W (1 GW)
    S = 10.0  # Photovoltaic cell area in m^2
    T = 12 * 3600  # Orbital period in s (12 hours)
    M = 7.35e22  # Lunar mass in kg
    R = 1738e3  # Lunar radius in m (1738 km)
    G = 6.67e-11  # Gravitational constant in kg^-1 m^3 s^-2

    # Step 1: Calculate the semi-major axis 'a' of the orbit using Kepler's Third Law.
    # a^3 = (G * M * T^2) / (4 * pi^2)
    a_cubed = (G * M * T**2) / (4 * math.pi**2)
    a = a_cubed**(1/3.0)

    # Step 2 & 3: Assume the total optical path length L is the major axis of the orbit, 2a.
    L = 2 * a

    # Step 4: The effective area of the horizontal cell for a vertical incident ray is S.
    S_eff = S

    # Step 5: Calculate the final power P' using the formula P' = P * S / (16 * pi * a^2)
    P_prime = (P * S_eff) / (4 * math.pi * L**2)

    # Step 6: Convert the result to microwatts.
    P_prime_microwatts = P_prime * 1e6

    # --- Output the calculation steps and results ---
    print("--- Calculation Steps ---")
    print(f"Given source power P = {P:.1e} W")
    print(f"Given cell area S = {S:.1f} m^2")
    print(f"Given orbital period T = {T:.0f} s")
    print(f"Given lunar mass M = {M:.2e} kg")
    print(f"Given gravitational constant G = {G:.2e} kg^-1 m^3 s^-2")
    
    print("\n1. Calculate semi-major axis 'a':")
    print(f"a = ((G * M * T^2) / (4 * pi^2))^(1/3)")
    print(f"a = (({G:.2e} * {M:.2e} * {T:.0f}^2) / (4 * {math.pi:.4f}^2))^(1/3)")
    print(f"a = {a:.4e} m")

    print("\n2. Determine total optical path length L:")
    print("Based on the problem's symmetries and implicit constraints, we assume L = 2a.")
    print(f"L = 2 * {a:.4e} m = {L:.4e} m")

    print("\n3. Calculate the final power P' on the cell:")
    print(f"P' = P * S / (4 * pi * L^2)")
    print(f"P' = ({P:.1e} * {S:.1f}) / (4 * {math.pi:.4f} * ({L:.4e})^2)")
    print(f"P' = {P_prime:.4e} W")

    print("\n4. Convert power to microwatts:")
    print(f"P' (uW) = {P_prime:.4e} W * 1,000,000 uW/W = {P_prime_microwatts:.4f} uW")

    print("\n--- Final Answer ---")
    print(f"The total power incident on the cell is {P_prime_microwatts:.1f} microwatts.")

calculate_incident_power()