import math

def solve_pandora_measurements():
    """
    Analyzes measurements from the star Pandora to find the likely error.
    
    Given:
    - Wavelength (λ): 400 nm
    - Spectral Radiance (B): 1.2e15 W/m²/sr/m
    - Surface Temperature (T): 9000 K

    Based on Planck's Law, the measured spectral radiance is too high for the
    estimated temperature. This script calculates the corrected temperature,
    assuming the wavelength and spectral radiance measurements are accurate.
    """
    
    # Scientific constants
    h = 6.62607015e-34  # Planck's constant in J·s
    c = 299792458       # Speed of light in m/s
    k_B = 1.380649e-23   # Boltzmann constant in J/K

    # Given measurements
    lambda_val = 400e-9  # Wavelength in meters
    B_measured = 1.2e15  # Measured spectral radiance in W/(m²·sr·m)
    
    # We solve Planck's Law for Temperature (T). The equation is:
    # T = (h*c) / (λ * k_B * ln( (2*h*c² / (λ⁵ * B)) + 1 ))

    try:
        # Calculate the argument for the natural logarithm
        log_argument = (2 * h * c**2) / (lambda_val**5 * B_measured) + 1
        
        if log_argument <= 0:
            print("Error: Invalid physical values, cannot compute logarithm.")
            return

        # Calculate the exponent term in Planck's law
        exponent_term = math.log(log_argument)
        
        # Calculate the corrected temperature
        T_correct = (h * c) / (lambda_val * k_B * exponent_term)
        
        # Round the result to the nearest integer as requested
        T_correct_rounded = int(round(T_correct))

        # The final answer format is Xv, where X is the quantity with the error
        # and v is its expected correct value.
        # X is Temperature (T)
        # v is the corrected value T_correct_rounded
        
        print("Final Equation solution:")
        print(f"Corrected Temperature (v) = {T_correct_rounded} K")
        print("The quantity most likely in error is Temperature (X) = T")
        print("\nFinal Answer (form: Xv):")
        print(f"T{T_correct_rounded}")

    except Exception as e:
        print(f"An error occurred during calculation: {e}")

# Execute the function to find the answer
solve_pandora_measurements()