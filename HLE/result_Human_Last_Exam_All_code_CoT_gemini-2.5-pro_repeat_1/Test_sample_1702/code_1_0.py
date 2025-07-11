import math

def calculate_temperature():
    """
    Calculates the temperature of the star Pandora based on its spectral radiance,
    using an approximation for Planck's Law suitable for a system without
    floating-point or advanced math functions.
    """
    # Constants
    h = 6.626e-34  # Planck's constant (J*s)
    c = 2.998e8    # Speed of light (m/s)
    k = 1.380e-23  # Boltzmann constant (J/K)
    
    # Given observations
    # Wavelength in meters (500 nm -> 500e-9 m)
    lambda_val = 500e-9
    # Spectral radiance (W/m^2*sr*m)
    B = 9.9e16

    # Step 1: Calculate the term 'y' from the rearranged Planck's Law.
    # y = (2 * h * c^2) / (B * λ^5)
    numerator_y = 2 * h * c**2
    denominator_y = B * lambda_val**5
    y = numerator_y / denominator_y

    # Step 2: Approximate x = ln(1+y) using the Taylor series: x ≈ y - y^2 / 2
    # This is valid because y is small.
    x = y - (y**2) / 2

    # Step 3: Calculate the constant factor from the exponent term.
    # C = hc / (λk)
    exponent_constant = (h * c) / (lambda_val * k)

    # Step 4: Calculate the final temperature T.
    # T = C / x
    temperature = exponent_constant / x

    # Step 5: Round the temperature to the nearest thousand Kelvin.
    # For example, 51723 K -> (51723 + 500) / 1000 = 52.
    rounded_temp_in_kK = int((temperature + 500) / 1000)

    # Output the final equation with the calculated values as required.
    print(f"The temperature T is calculated from the formula: T = (hc / (λk)) / x")
    print(f"where x is an approximation for ln(1 + (2hc^2)/(Bλ^5)).")
    print(f"Calculated value for the constant (hc / (λk)): {exponent_constant:.0f} K")
    print(f"Calculated value for x: {x:.6f}")
    print(f"Final equation: T = {exponent_constant:.0f} / {x:.6f}")
    print(f"Resulting temperature T = {temperature:.0f} K")
    print(f"\nThe temperature of Pandora, rounded to the nearest thousand Kelvin, is: {rounded_temp_in_kK}")

calculate_temperature()