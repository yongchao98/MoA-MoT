import sympy

def solve_mass_ratio():
    """
    Calculates the leading-order asymptotic mass ratio for the CP(N-1) model.

    In the large-N limit, the mass spectrum of the CP(N-1) model includes solitonic
    excitations with mass M_k = (N / (2 * pi)) * m * sin(k * pi / N), where m is the
    mass of the fundamental mesons.

    This function calculates the ratio of the mass of the lightest soliton (k=1)
    to the mass of the subsequent higher excitation (a meson of mass m) as N
    approaches infinity.
    """
    # Define the symbolic variables. N is the model parameter, and m is the
    # mass of the fundamental meson, which serves as our mass scale.
    N, m = sympy.symbols('N m', positive=True)

    # The lightest solitonic excitation has topological charge k=1.
    k = 1
    mass_lightest_soliton = (N / (2 * sympy.pi)) * m * sympy.sin(k * sympy.pi / N)

    # The subsequent higher excitation in the spectrum has mass m.
    mass_next_excitation = m

    # Form the ratio of the two masses. The mass 'm' cancels out.
    mass_ratio_expression = mass_lightest_soliton / mass_next_excitation

    # Calculate the asymptotic limit of this ratio as N approaches infinity.
    asymptotic_ratio = sympy.limit(mass_ratio_expression, N, sympy.oo)

    # The final result is a rational number. We extract the numerator and
    # denominator to explicitly show the numbers in the final equation,
    # satisfying the problem requirements.
    # Final Equation: Ratio = Numerator / Denominator
    if asymptotic_ratio.is_Rational:
        numerator = asymptotic_ratio.p
        denominator = asymptotic_ratio.q
    else:
        # Fallback for unexpected result types
        numerator, denominator = sympy.fraction(asymptotic_ratio)

    print(f"The formula for the mass ratio is: {mass_ratio_expression}")
    print(f"The calculated asymptotic limit as N -> infinity is: {asymptotic_ratio}")
    print("\nThe final equation for the ratio is:")
    print(f"Ratio = {numerator} / {denominator}")
    print("\nComponents of the final equation:")
    print(f"Numerator: {numerator}")
    print(f"Denominator: {denominator}")

if __name__ == '__main__':
    solve_mass_ratio()