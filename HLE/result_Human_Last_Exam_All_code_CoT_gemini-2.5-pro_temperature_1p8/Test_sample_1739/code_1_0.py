from fractions import Fraction

def calculate_frequency_correction_terms():
    """
    Calculates and prints the terms of the nonlinear frequency correction coefficient
    for the Rayleigh-Plesset equation.

    The dimensionless frequency correction is of the form:
    omega / omega_0 = 1 - C(gamma) * a_R^2
    where C(gamma) is a polynomial in the polytropic index gamma.
    """

    # The unsimplified form of the coefficient C(gamma) derived from a
    # two-time scale analysis is:
    # C(gamma) = (72 * gamma^2 + 45 * gamma + 6) / 48
    numerator_coeffs = [72, 45, 6]
    denominator = 48

    # Calculate the simplified fractional coefficients for C(gamma) = C2*gamma^2 + C1*gamma + C0
    C2 = Fraction(numerator_coeffs[0], denominator)
    C1 = Fraction(numerator_coeffs[1], denominator)
    C0 = Fraction(numerator_coeffs[2], denominator)

    print("The nonlinear correction to the linear oscillation frequency is expressed as:")
    print("omega/omega_0 = 1 - C(gamma) * a_R^2")
    print("where a_R is the amplitude of the radius oscillation and C(gamma) is a coefficient that depends on the polytropic index gamma.")
    print("\nThe derived polynomial for the coefficient C(gamma) is:")
    
    # Final equation showing all the calculated numbers
    print(f"C(gamma) = ({C2.numerator}/{C2.denominator})*gamma^2 + ({C1.numerator}/{C1.denominator})*gamma + ({C0.numerator}/{C0.denominator})")

    print("\nThe nonlinear correction consists of three terms:")
    print(f"1st term: ({C2.numerator}/{C2.denominator})*gamma^2")
    print(f"2nd term: ({C1.numerator}/{C1.denominator})*gamma")
    print(f"3rd term: {C0.numerator}/{C0.denominator}")

    # The user asked to calculate the 3rd term of the correction.
    third_term_value = C0
    print(f"\nThe value of the 3rd term is: {third_term_value}")

# Execute the calculation
calculate_frequency_correction_terms()