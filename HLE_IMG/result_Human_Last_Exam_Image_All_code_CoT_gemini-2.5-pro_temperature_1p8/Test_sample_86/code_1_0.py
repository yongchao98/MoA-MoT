import sympy as sp

def print_energy_spectrum():
    """
    This function prints the derived energy spectrum for the perturbed harmonic oscillator.
    The final equation is constructed and printed step-by-step to show all the numerical coefficients.
    """
    # Define symbols for mathematical representation, though they won't be used for computation here.
    n, hbar, u, m, w0 = sp.symbols('n hbar u m omega_0')

    # The problem requires printing the numbers in the final equation.
    # We will construct the string for the final equation piece by piece.
    
    term1 = "(n + 1/2) * hbar * omega_0"
    
    # The correction term is (n + 1/2) * u * hbar^2 / (8 * m^2 * omega_0^2)
    # The coefficients are 1 for the numerator and 8 for the denominator.
    numerator_constant = 1
    denominator_constant = 8
    
    term2_factor = "(n + 1/2)"
    term2_fraction = f"({numerator_constant} * u * hbar**2) / ({denominator_constant} * m**2 * omega_0**2)"
    
    term2 = f"{term2_factor} * {term2_fraction}"
    
    final_equation = f"E_n = {term1} + {term2}"

    print("The self-energy diagram predicts the following energy spectrum, up to an overall constant E_0:")
    print("The final equation is:")
    print(final_equation)
    
    print("\nHere are the explicit numerical constants from the formula above:")
    print("Constant term added to n: 1/2")
    print("Numerator constant in the correction term fraction: 1")
    print("Denominator constant in the correction term fraction: 8")

if __name__ == '__main__':
    print_energy_spectrum()