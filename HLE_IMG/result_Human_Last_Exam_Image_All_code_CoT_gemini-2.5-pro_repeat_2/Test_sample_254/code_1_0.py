import math

def solve_power_on_cell():
    """
    Calculates the power incident on a photovoltaic cell on the Moon's surface
    from a light source, relayed by two satellites.
    """
    # --- Given Data ---
    # Power of the isotropic source in Watts (1 GW)
    P = 1.0e9
    # Mass of the Moon in kg
    M = 7.35e22
    # Radius of the Moon in meters (1738 km)
    R = 1738e3
    # Gravitational constant in kg^-1 m^3 s^-2
    G = 6.67e-11
    # Orbital period in seconds (12 hours)
    T = 12 * 3600
    # Area of the photovoltaic cell in m^2
    S = 10.0

    print("Step 1: Calculate the orbital radius from Kepler's Third Law.")
    # For a circular orbit, the semi-major axis 'a' is the radius 'r'.
    # T^2 = (4 * pi^2 / (G * M)) * r^3  =>  r = (G * M * T^2 / (4 * pi^2))^(1/3)
    r_cubed = (G * M * T**2) / (4 * math.pi**2)
    r = r_cubed**(1/3)
    print(f"The orbital radius is r = {r:.2e} m or {r/1000:.1f} km.\n")

    print("Step 2: Calculate the satellite's altitude.")
    # The altitude 'h' is the orbital radius minus the Moon's radius.
    # This altitude corresponds to the length of the vertical light paths AX and YB.
    h = r - R
    d_YB = h
    print(f"The satellite altitude is h = r - R = {h:.2e} m or {h/1000:.1f} km.")
    print(f"Due to the symmetric setup, the distance d_YB = h.\n")

    print("Step 3: Calculate the power P' on the cell.")
    # The intensity at the final destination B is I_B = P / (2 * pi * d_YB^2).
    # The power on the cell is P' = I_B * S.
    # P' = (P / (2 * pi * d_YB^2)) * S
    intensity_at_B = P / (2 * math.pi * d_YB**2)
    power_incident_cell = intensity_at_B * S
    
    # Convert to microwatts
    power_in_microwatts = power_incident_cell * 1e6
    
    print("The final power is calculated using the formula: P' = (P / (2 * pi * d_YB^2)) * S")
    print("Substituting the values:")
    print(f"P' = ({P:.1e} W / (2 * pi * ({d_YB:.2e} m)^2)) * {S:.1f} m^2")
    print(f"P' = {power_incident_cell:.6e} W")
    print(f"\nFinal Answer: The total power incident on the cell is {power_in_microwatts:.1f} microwatts.")
    
    return power_in_microwatts

# Run the calculation and print the final result
final_power = solve_power_on_cell()
# The final answer in the requested format
# print(f"\n<<<{final_power:.1f}>>>")