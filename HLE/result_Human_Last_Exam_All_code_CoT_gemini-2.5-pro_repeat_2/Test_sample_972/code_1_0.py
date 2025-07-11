import sympy

def solve_wave_amplitude():
    """
    This function calculates and displays the symbolic formula for the electric field
    amplitude of a wave exiting a time-varying slab.

    The derivation is based on the following steps:
    1.  The medium is impedance-matched to vacuum (Z=Z_0), so there are no reflections.
    2.  The adiabatic invariant for wave propagation in a slowly varying medium,
        S/omega = constant, is applied. Since Z is constant, this implies E^2/omega = constant.
    3.  The change in frequency (omega) is calculated by finding the transit time through
        the slab. The wave speed is v(t) = c/n(t) where n(t) = alpha*t + beta.
    4.  The integration of dx/dt = v(t) gives the relation between input and output frequencies:
        omega_out = omega_in * exp(-alpha * L / c).
    5.  Combining these facts gives the final amplitude E_L.
    """

    # Define the symbols for the parameters in the problem.
    # A: Amplitude of the incident wave
    # alpha: Rate of change of the material property
    # L: Length of the slab
    # c: Speed of light in vacuum (a constant, but represented symbolically)
    A, alpha, L, c = sympy.symbols('A alpha L c')

    # The final equation for the amplitude E_L at x=L is:
    # E_L = A * sqrt(omega_out / omega_in)
    # E_L = A * sqrt(exp(-alpha * L / c))
    # E_L = A * exp(-alpha * L / (2 * c))

    # We construct the final expression with explicit numbers as requested by the prompt.
    # The numbers in the final equation's exponent are -1 and 2.
    numerator_coeff = -1
    denominator_coeff = 2

    # Build the final equation string showing the numbers explicitly.
    final_equation_str = f"E_L = {A} * exp(({numerator_coeff} * {alpha} * {L}) / ({denominator_coeff} * {c}))"

    print("The amplitude of the electric field at the rightmost boundary of the slab is given by the formula:")
    print(final_equation_str)

    # For better mathematical representation, we can also use sympy's pretty printing.
    final_expression_sympy = A * sympy.exp(numerator_coeff * alpha * L / (denominator_coeff * c))
    print("\nIn standard mathematical notation:")
    sympy.pprint(final_expression_sympy, use_unicode=True)


if __name__ == '__main__':
    solve_wave_amplitude()
