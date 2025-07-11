import sys

def get_function_name_and_formula():
    """
    This function provides the name and a common formula for the space-time,
    double Fourier transform of the generalized pair correlation function
    in nuclear criticality.
    """
    
    # 1. Identify and state the common name of the function.
    function_name = "Power Spectral Density (PSD) of the Neutron Noise"
    
    print("The space-time, double Fourier transform of the generalized pair correlation function is commonly called the:")
    print(f"'{function_name}'\n")

    # 2. Define and display the formula.
    # This formula represents the PSD in a simple one-group diffusion model.
    numerator = "C"
    denominator = "((alpha + D*v*k**2)**2 + omega**2)"
    formula_str = f"S(k, omega) = {numerator} / {denominator}"
    
    print("A common theoretical form for this function is:")
    print(formula_str)
    print("\n-------------------------------------------")

    # 3. Explain the terms in the formula.
    print("Where the terms are defined as:")
    variables = {
        "S(k, omega)": "The Power Spectral Density as a function of wavenumber (k) and angular frequency (omega).",
        "C": "A constant proportional to the fission rate and the Diven factor (X2 * F / Lambda**2).",
        "alpha": "The prompt neutron decay constant, also known as Rossi-alpha.",
        "D": "The neutron diffusion coefficient.",
        "v": "The average neutron speed.",
        "k": "The wavenumber, representing the spatial frequency.",
        "omega": "The angular frequency, representing the temporal frequency.",
        "**2": "This indicates a variable is raised to the power of 2 (squared)."
    }
    
    for term, description in variables.items():
        print(f"- {term}: {description}")

    print("-------------------------------------------\n")

    # 4. Fulfill the request to output the numbers in the equation.
    # The only numbers explicitly written in this formula are the exponents.
    print("The explicit numbers appearing in the formula are the exponents:")
    print("The number 2 appears as an exponent on:")
    print("- The term k (wavenumber)")
    print("- The term omega (angular frequency)")
    print("- The entire term (alpha + D*v*k**2)")


if __name__ == '__main__':
    get_function_name_and_formula()
