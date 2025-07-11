import fractions

def calculate_relativistic_shift():
    """
    Calculates the second-order energy shift due to the relativistic
    kinetic energy correction for a hydrogen atom.
    """
    # Define the quantum numbers for the given state
    n = 3
    l = 2

    print(f"Calculating the shift for n={n}, l={l}.")
    print("-" * 30)

    # The calculation will be done using the fractions module for exact results.
    n_f = fractions.Fraction(n)
    l_f = fractions.Fraction(l)

    # The formula for the shift is:
    # ΔE = - (m * c^2 * α^6) / (4 * n^3) * [ 1/(l+0.5)^2 + 4/(l+0.5) - 1 - 2/(n*(l+0.5)) ]
    # We will calculate the numerical coefficient for this expression.

    # Calculate l + 1/2
    l_plus_half = l_f + fractions.Fraction(1, 2)

    # Calculate the term inside the square brackets
    term1 = 1 / (l_plus_half**2)
    term2 = 4 / l_plus_half
    term3 = fractions.Fraction(-1)
    term4 = -2 / (n_f * l_plus_half)
    bracket_term = term1 + term2 + term3 + term4

    # Calculate the full numerical coefficient
    coefficient = -fractions.Fraction(1, 4 * n**3) * bracket_term

    # Extract the numerator and denominator of the final coefficient
    final_numerator = coefficient.numerator
    final_denominator = coefficient.denominator

    # The final energy shift is ΔE = (numerator/denominator) * m * c^2 * α^6
    # As requested, we will output each number in the final equation.
    print("The second-order energy shift due to the relativistic kinetic energy correction is:")
    print(f"ΔE = ({final_numerator} / {final_denominator}) * m * c^2 * α^6")
    
    print("\nTo be explicit, the numbers in the final equation are:")
    print(f"  - Numerator of the coefficient: {final_numerator}")
    print(f"  - Denominator of the coefficient: {final_denominator}")
    print(f"  - The fundamental constants are electron mass 'm' and the speed of light 'c'.")
    print(f"  - The power of c is 2.")
    print(f"  - The fine-structure constant is 'α'.")
    print(f"  - The power of α is 6.")

if __name__ == '__main__':
    calculate_relativistic_shift()
