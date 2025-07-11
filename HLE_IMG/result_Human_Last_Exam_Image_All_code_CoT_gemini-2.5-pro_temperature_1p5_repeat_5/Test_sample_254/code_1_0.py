import math

def calculate_power_on_cell():
    """
    Calculates the power incident on a photovoltaic cell on the Moon's surface
    from a light source reflected by two satellites.
    """
    # --- Given Data (in SI units) ---
    # Luminosity of the isotropic point source (Watts)
    P = 1.0 * 10**9  # 1 GW
    # Mass of the Moon (kg)
    M = 7.35 * 10**22
    # Radius of the Moon (meters)
    R = 1738 * 1000
    # Universal Gravitational Constant (m^3 kg^-1 s^-2)
    G = 6.67 * 10**-11
    # Orbital period of the satellites (seconds)
    T = 12 * 3600
    # Area of the photovoltaic cell (m^2)
    S = 10.0

    # --- Step 1: Calculate the semi-major axis 'a' of the orbit ---
    # From Kepler's Third Law: a = (G * M * T^2 / (4 * pi^2))^(1/3)
    a = ((G * M * T**2) / (4 * math.pi**2))**(1/3)

    # --- Step 2: Calculate the altitude 'h' of the satellite ---
    # Assuming a circular orbit for the specified reflection geometry to be possible.
    # The altitude 'h' is the distance from the source A to satellite X.
    h = a - R

    # --- Step 3: Calculate the power P' incident on the cell ---
    # The power is given by the formula P' = (P * S) / (4 * pi * h^2)
    P_prime_watts = (P * S) / (4 * math.pi * h**2)

    # Convert the result to microwatts
    P_prime_microwatts = P_prime_watts * 10**6

    # --- Print the detailed calculation steps ---
    print("Problem solving steps:")
    print("1. Calculate the semi-major axis (a) of the orbit using Kepler's Third Law.")
    print(f"   a = (G * M * T^2 / (4 * pi^2))^(1/3)")
    print(f"   a = (({G:.2e} kg^-1 m^3 s^-2 * {M:.2e} kg * ({T} s)^2) / (4 * {math.pi:.4f}^2))^(1/3)")
    print(f"   a = {a:.2f} m")
    print("\n2. Calculate the satellite's altitude (h) assuming a circular orbit.")
    print(f"   h = a - R = {a:.2f} m - {R:.2f} m = {h:.2f} m")
    print("\n3. Calculate the power (P') on the cell.")
    print(f"   P' = (P * S) / (4 * pi * h^2)")
    print(f"   P' = ({P:.1e} W * {S} m^2) / (4 * {math.pi:.4f} * ({h:.2f} m)^2)")
    print(f"   P' = {P_prime_watts:.6e} W")
    print("\n4. Convert power to microwatts and round to one decimal place.")
    print(f"   P' (µW) = {P_prime_watts:.6e} W * 1,000,000 = {P_prime_microwatts:.1f} µW")


calculate_power_on_cell()
<<<41.0>>>