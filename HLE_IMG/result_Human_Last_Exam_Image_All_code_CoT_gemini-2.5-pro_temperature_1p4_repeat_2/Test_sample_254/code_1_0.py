import math

def calculate_incident_power():
    """
    Calculates the total power incident on the photovoltaic cell at point B.
    """
    # --- Given Data ---
    # Power of the isotropic source in Watts
    P = 1e9
    # Area of the photovoltaic cell at B in m^2
    S = 10
    # Orbital period of the satellites in seconds
    T = 12 * 3600

    # --- Physical Constants ---
    # Gravitational constant in m^3 kg^-1 s^-2
    G = 6.67e-11
    # Mass of the Moon in kg
    M = 7.35e22

    # Step 1: Calculate the semi-major axis 'a' of the orbit using Kepler's Third Law.
    # T^2 = (4 * pi^2 * a^3) / (G * M)
    a_cubed = (G * M * T**2) / (4 * math.pi**2)
    a = a_cubed**(1/3)

    # Step 2: Determine the total light path length.
    # For the specific geometry described, the total path length d_total = 4 * a.
    d_total = 4 * a

    # Step 3: Determine the effective area of the cell.
    # The geometry ensures the final ray is normal to the horizontal cell, so the incidence angle is 0.
    # The effective area is equal to the cell's actual area S.
    S_eff = S

    # Step 4: Calculate the final power P' incident on the cell.
    # The power from an isotropic source spreads over a sphere of area 4 * pi * d^2.
    # P' = P_source * Area_cell / (4 * pi * total_distance^2)
    P_prime = (P * S_eff) / (4 * math.pi * d_total**2)
    
    # Step 5: Convert the answer from Watts to microwatts.
    P_prime_microwatts = P_prime * 1e6
    
    print("--- Calculation Steps ---")
    print(f"1. The semi-major axis of the orbit 'a' is {a/1000:.1f} km.")
    print(f"2. The total path length of light d_total = 4 * a = {d_total/1000:.1f} km.")
    print(f"3. The final power P' is calculated using the formula:")
    print(f"   P' = P * S / (4 * pi * d_total^2)")
    print(f"   P' = {P:.0e} W * {S} m^2 / (4 * pi * ({d_total:.2e} m)^2)")
    print(f"   P' = {P_prime:.3e} W")
    print(f"4. Converting to microwatts: {P_prime:.3e} W * 1,000,000 uW/W = {P_prime_microwatts:.1f} uW.")
    print("------------------------")
    print(f"The total power incident on the cell is {P_prime_microwatts:.1f} microwatts.")

calculate_incident_power()
<<<1.3>>>