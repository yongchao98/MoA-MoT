import fractions

def calculate_relativistic_shift():
    """
    Calculates the relativistic energy shift for a hydrogen atom state.
    
    The energy shift is given by the first-order perturbation theory formula:
    ΔE = <n,l,m| H' |n,l,m> = -(α^4*m*c^2 / (2*n^4)) * (n/(l+1/2) - 3/4)
    """
    
    # Given principal and angular momentum quantum numbers
    n = 3
    l = 2
    
    print(f"Calculating the relativistic energy shift for n={n}, l={l}.")

    # We use the `fractions` module to maintain precision in our calculations.
    n_f = fractions.Fraction(n)
    l_f = fractions.Fraction(l)

    # Calculate the numerical coefficient of the expression: α^4 * m * c^2
    # The formula for the coefficient is: [-1 / (2 * n^4)] * [n / (l + 1/2) - 3/4]
    
    # Part 1 of the coefficient: [-1 / (2 * n^4)]
    coeff_part1 = -fractions.Fraction(1, 2 * n**4)
    
    # Part 2 of the coefficient: [n / (l + 1/2) - 3/4]
    coeff_part2 = n_f / (l_f + fractions.Fraction(1, 2)) - fractions.Fraction(3, 4)
    
    # Total coefficient
    total_coeff = coeff_part1 * coeff_part2
    
    numerator = total_coeff.numerator
    denominator = total_coeff.denominator

    # Display the final derived equation for the energy shift
    print("\nThe final expression for the energy shift is:")
    print(f"ΔE = ({numerator}/{denominator}) * α^4 * m * c^2")
    
    # Outputting the numbers in the final equation per the user's request
    print("\nThe numbers in the final equation are:")
    print(f"Numerator of the coefficient: {numerator}")
    print(f"Denominator of the coefficient: {denominator}")
    print("Exponent of α (fine-structure constant): 4")
    print("Exponent of m (electron mass): 1")
    print("Exponent of c (speed of light): 2")

calculate_relativistic_shift()