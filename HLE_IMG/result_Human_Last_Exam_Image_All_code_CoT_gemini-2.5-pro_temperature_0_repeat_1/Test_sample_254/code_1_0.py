import math

def solve_power_on_cell():
    """
    Calculates the total power incident on a photovoltaic cell on the Moon's surface
    from a light source via two satellite reflections.
    """
    # Given data
    P_GW = 1  # GW
    T_hours = 12  # hours
    S_m2 = 10  # m^2
    M_kg = 7.35 * 10**22  # kg
    R_km = 1738  # km
    G = 6.67 * 10**-11  # kg^-1 m^3 s^-2

    # Convert to SI units
    P = P_GW * 10**9  # Watts
    T = T_hours * 3600  # seconds
    S = S_m2 # m^2
    R = R_km * 10**3  # meters
    pi = math.pi

    # Step 1: Calculate the semi-major axis 'a' of the orbit
    # From Kepler's Third Law: a = (G * M * T^2 / (4 * pi^2))^(1/3)
    a_cubed = (G * M * T**2) / (4 * pi**2)
    a = a_cubed**(1/3)

    # Step 2: Determine the distance from the source to the first satellite (d_AX)
    # Based on the problem's constraints, we assume a special configuration
    # where the satellite's orbital radius r_X equals the semi-major axis a.
    # The satellite X is at the zenith of point A on the surface.
    r_X = a
    d_AX = r_X - R

    # Step 3: Calculate the final power P' on the cell
    # The power formula simplifies under the assumption of a symmetric reflection path to:
    # P' = (P * S) / (4 * pi * d_AX^2)
    P_prime = (P * S) / (4 * pi * d_AX**2)

    # Convert the final answer to microwatts
    P_prime_microwatts = P_prime * 10**6

    # Print the calculation steps and the final result
    print("Calculation Steps:")
    print(f"1. The semi-major axis of the orbit is calculated using Kepler's Third Law.")
    print(f"   a = (G * M * T^2 / (4 * pi^2))^(1/3)")
    print(f"   a = ({G:.2e} * {M:.2e} * {T}^2 / (4 * {pi:.4f}^2))^(1/3) = {a:.4e} m")
    print("\n2. The distance from the source A to satellite X is determined.")
    print(f"   Assuming the satellite is at a distance from the Moon's center equal to the semi-major axis, r_X = a.")
    print(f"   d_AX = r_X - R = {a:.4e} m - {R:.4e} m = {d_AX:.4e} m")
    print("\n3. The power on the cell is calculated.")
    print(f"   P' = (P * S) / (4 * pi * d_AX^2)")
    print(f"   P' = ({P:.2e} W * {S} m^2) / (4 * {pi:.4f} * ({d_AX:.4e} m)^2) = {P_prime:.4e} W")
    print("\n4. The final power is converted to microwatts.")
    print(f"   P' = {P_prime:.4e} W * 10^6 uW/W = {P_prime_microwatts:.1f} uW")

solve_power_on_cell()
<<<41.0>>>