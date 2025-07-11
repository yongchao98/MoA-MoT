import math

def calculate_power_on_cell():
    """
    Calculates the power incident on a photovoltaic cell on the Moon's surface
    after light from a source is reflected by two satellites.
    """

    # --- Given Data (in SI units) ---
    # Source Power (P) in Watts
    P = 1 * 10**9
    # Orbital Period (T) in seconds
    T = 12 * 3600
    # Lunar Mass (M) in kg
    M = 7.35 * 10**22
    # Lunar Radius (R) in meters
    R = 1738 * 1000
    # Gravitational Constant (G) in m^3 kg^-1 s^-2
    G = 6.67 * 10**-11
    # Area of the photovoltaic cell (S) in m^2
    S = 10
    
    # --- Step 1: Calculate the orbital radius (semi-major axis 'a') ---
    # Using Kepler's Third Law: a^3 = (G * M * T^2) / (4 * pi^2)
    a_cubed = (G * M * T**2) / (4 * math.pi**2)
    a = a_cubed**(1/3)
    
    # --- Step 2: Calculate the square of the total optical path length (L^2) ---
    # Based on the geometry of two reflections creating a double virtual source,
    # the squared distance from the final virtual source A'' to B is L^2 = d(A'',B)^2.
    # For the symmetric 90-degree configuration, this simplifies to:
    # L_sq = (2*a - R)^2 + R^2
    L_sq = (2 * a - R)**2 + R**2

    # --- Step 3: Calculate the final power P' incident on the cell ---
    # P' = P_source * Area_cell / (4 * pi * L_sq)
    P_prime_watts = (P * S) / (4 * math.pi * L_sq)
    
    # Convert power to microwatts (1 W = 1,000,000 uW)
    P_prime_microwatts = P_prime_watts * 1_000_000

    # --- Output the results ---
    print("--- Calculation Steps ---")
    print(f"1. Calculated orbital radius (a): {a:.2f} meters")
    print("\n2. The final power P' is calculated using the formula: P' = (P * S) / (4 * pi * L^2)")
    print("   Where:")
    print(f"   - P (source power) = {P:.1e} W")
    print(f"   - S (cell area) = {S} m^2")
    print(f"   - pi = {math.pi:.4f}")
    print(f"   - L^2 (squared optical path length) = {L_sq:.4e} m^2")
    
    print("\n3. Plugging in the values:")
    print(f"   P' = ({P:.1e} * {S}) / (4 * {math.pi:.4f} * {L_sq:.4e})")
    print(f"   P' = {P_prime_watts:.4e} Watts")

    print("\n--- Final Answer ---")
    print(f"The total power P' incident on the cell is {P_prime_microwatts:.1f} microwatts.")

calculate_power_on_cell()
<<<7.0>>>