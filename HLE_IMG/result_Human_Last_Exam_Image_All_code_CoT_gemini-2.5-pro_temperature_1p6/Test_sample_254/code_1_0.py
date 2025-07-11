import math

def solve_power_on_cell():
    """
    Calculates the power incident on a photovoltaic cell on the Moon's surface
    after light from a source is reflected by two satellites.
    """
    # Given data in SI units
    P = 1.0e9  # Luminosity of the source in Watts
    M = 7.35e22  # Mass of the Moon in kg
    R = 1738e3   # Radius of the Moon in meters
    G = 6.67e-11 # Gravitational constant in m^3 kg^-1 s^-2
    T = 12 * 3600  # Orbital period in seconds
    S_cell = 10.0  # Area of the photovoltaic cell in m^2

    print("--- Given Data ---")
    print(f"Source Luminosity (P): {P:.1e} W")
    print(f"Lunar Mass (M): {M:.2e} kg")
    print(f"Lunar Radius (R): {R:.0f} m")
    print(f"Orbital Period (T): {T:.0f} s")
    print(f"Cell Area (S): {S_cell:.1f} m^2")
    print("-" * 20 + "\n")

    # Step 1: Calculate the semi-major axis 'a' of the orbit using Kepler's Third Law
    # T^2 = (4 * pi^2 / (G * M)) * a^3
    a_cubed = (G * M * T**2) / (4 * math.pi**2)
    a = a_cubed**(1/3)
    
    print("--- Step 1: Calculate Orbit's Semi-Major Axis 'a' ---")
    print(f"Equation: a = (G * M * T^2 / (4 * pi^2))^(1/3)")
    print(f"a = ({G:.2e} * {M:.2e} * {T**2:.2e} / (4 * {math.pi**2:.2f}))^(1/3)")
    print(f"Result: a = {a:.1f} m\n")

    # Step 2: Determine the pericenter distance r_p using the geometric constraint
    # (r_p - R) * (r_a - R) = R^2  and  r_p + r_a = 2a
    # This leads to a quadratic equation for r_p: r_p^2 - 2*a*r_p + 2*a*R = 0
    # Solving for r_p: r_p = a - sqrt(a^2 - 2*a*R) (we take the smaller root for pericenter)
    discriminant = a**2 - 2 * a * R
    if discriminant < 0:
        print("Error: Orbit not possible with these parameters (negative discriminant).")
        return

    r_p = a - math.sqrt(discriminant)
    
    print("--- Step 2: Calculate Pericenter Distance 'r_p' ---")
    print(f"From the geometric constraint (h_p * h_a = R^2), we solve:")
    print(f"Equation: r_p^2 - 2*a*r_p + 2*a*R = 0")
    print(f"r_p = {a:.1f} - sqrt({a**2:.2e} - 2 * {a:.1f} * {R:.0f})")
    print(f"Result: r_p = {r_p:.1f} m\n")
    
    if r_p <= R:
        print("Error: Satellite would crash into the Moon.")
        return

    # Step 3: Calculate the altitude 'h' of satellite X at the pericenter
    # This altitude is the distance d_AX from the source A to satellite X.
    h = r_p - R
    
    print("--- Step 3: Calculate Altitude 'h' of Satellite X ---")
    print(f"The satellite X is at the pericenter, at the zenith of point A.")
    print(f"Equation: h = r_p - R")
    print(f"h = {r_p:.1f} - {R:.0f}")
    print(f"Result: h = {h:.1f} m\n")
    
    # Step 4: Calculate the power P' incident on the cell
    # The final simplified formula is P' = P * S_cell / (4 * pi * h^2)
    P_prime = (P * S_cell) / (4 * math.pi * h**2)

    print("--- Step 4: Calculate Final Power P' on the Cell ---")
    print("Due to reciprocity, P' = P * S_cell / (4 * pi * h^2)")
    print(f"P' = {P:.1e} * {S_cell:.1f} / (4 * {math.pi:.2f} * {h:.1f}^2)")
    print(f"Result: P' = {P_prime:.6f} W\n")

    # Step 5: Convert the result to microwatts and round to one decimal place
    P_prime_microwatts = P_prime * 1e6
    
    print("--- Step 5: Final Answer in Microwatts ---")
    print(f"Power in microwatts = {P_prime_microwatts:.1f} uW")
    
    # Final answer in specified format
    print(f"\n<<<{P_prime_microwatts:.1f}>>>")

solve_power_on_cell()