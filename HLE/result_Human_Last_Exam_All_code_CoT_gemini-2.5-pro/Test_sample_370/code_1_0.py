import sympy

def calculate_cross_section():
    """
    Calculates and prints the total cross section for psi-psi scattering
    in the high-energy limit.
    """

    # Define symbols for the equation
    g, E, pi, sigma = sympy.symbols('g E pi sigma')

    # Explain the result based on the derivation
    print("Based on the principles of Quantum Field Theory, the total cross section σ for the scattering")
    print("of two fermions (ψψ -> ψψ) via pseudoscalar interaction is calculated.")
    print("The calculation is performed in the high-energy limit (E >> m, M).\n")

    # --- Coefficients from the derivation ---
    
    # The integral of the spin-averaged matrix element squared over the
    # kinematically allowed solid angle for identical particles gives a factor of 3*g^4.
    numerator_coeff = 3

    # The flux and phase space factors in the cross-section formula contribute 1 / (32 * pi * s).
    denominator_coeff_s = 32
    
    # In the high-energy limit, s = 4 * E^2.
    s_factor = 4
    
    # The total numerical coefficient in the denominator is 32 * 4.
    denominator_coeff_total = denominator_coeff_s * s_factor

    # --- Printing the final equation ---
    
    print("The final result for the total cross section is derived as follows:")
    print(f"1. The integral over the phase space gives a numerator of {numerator_coeff} * g^4.")
    print(f"2. The denominator from the cross-section formula is ({denominator_coeff_s} * π * s).")
    print(f"3. In the high-energy limit, s is replaced by ({s_factor} * E^2).")
    print(f"4. The final denominator becomes {denominator_coeff_s} * {s_factor} * π * E^2 = {denominator_coeff_total} * π * E^2.\n")

    # Display the final equation using the calculated coefficients
    # This ensures we output each number as requested.
    final_equation = sympy.Eq(sigma, (numerator_coeff * g**4) / (denominator_coeff_total * pi * E**2))

    print("The total cross section σ is:")
    # Pretty print the equation
    sympy.pprint(final_equation, use_unicode=True)

if __name__ == '__main__':
    calculate_cross_section()