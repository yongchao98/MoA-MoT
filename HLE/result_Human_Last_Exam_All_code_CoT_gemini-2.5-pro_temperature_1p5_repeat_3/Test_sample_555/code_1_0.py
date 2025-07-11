from fractions import Fraction

def solve_dice_problem():
    """
    This function calculates the minimal expected number of rolls to generate a random
    digit (0-9) using a 7-sided die.

    It sets up and solves a system of linear equations representing the states of the
    optimal strategy. The states are defined by the number of "leftover" outcomes
    from a previous step.
    """

    # We have the following system of equations:
    # E = 2 + (9/49) * E9
    # E9 = 1 + (1/21) * E3
    # E3 = 1 + (1/21) * E1
    # E1 = 2 + (9/49) * E9

    # Let's use substitution to solve for E9, E3, E1.
    # From E1 = 2 + (9/49)E9, substitute into E3's equation.
    # E3 = 1 + (1/21) * (2 + (9/49)E9)
    # E3 = 1 + 2/21 + 9/(21*49)E9 = 23/21 + 3/(7*49)E9 = 23/21 + 3/343 * E9
    c_E1_E9_num, c_E1_E9_den = 9, 49
    
    c_E3_E1_num, c_E3_E1_den = 1, 21
    E3_const_num = 1 * c_E3_E1_den + c_E3_E1_num * 2
    E3_const_den = c_E3_E1_den
    c_E3_E9_num = c_E3_E1_num * c_E1_E9_num
    c_E3_E9_den = c_E3_E1_den * c_E1_E9_den
    E3_const = Fraction(E3_const_num, E3_const_den)
    c_E3_E9 = Fraction(c_E3_E9_num, c_E3_E9_den)

    # Now substitute the expression for E3 into E9's equation.
    # E9 = 1 + (1/21) * E3 = 1 + (1/21) * (23/21 + (3/343)E9)
    # E9 = 1 + 23/441 + 3/(21*343)E9 = 464/441 + 1/2401 * E9
    c_E9_E3_num, c_E9_E3_den = 1, 21
    E9_const_num = 1 * E3_const.denominator + c_E9_E3_num * E3_const.numerator
    E9_const_den = c_E9_E3_den * E3_const.denominator
    c_E9_E9_num = c_E9_E3_num * c_E3_E9.numerator
    c_E9_E9_den = c_E9_E3_den * c_E3_E9.denominator
    E9_const = Fraction(E9_const_num, E9_const_den)
    c_E9_E9 = Fraction(c_E9_E9_num, c_E9_E9_den)

    # Solve for E9: E9 * (1 - c_E9_E9) = E9_const
    E9 = E9_const / (1 - c_E9_E9)
    
    # Now that we have E9, we can find the total expected rolls E.
    # E = 2 + (9/49) * E9
    c_E_E9_num, c_E_E9_den = 9, 49
    c_E_E9 = Fraction(c_E_E9_num, c_E_E9_den)
    E_const = Fraction(2)
    E = E_const + c_E_E9 * E9

    print("The optimal strategy leads to a system of linear equations for the expected number of rolls.")
    print("Let E be the total expected rolls, and E_n be the expected additional rolls given n leftover outcomes.")
    
    # Print equation for E
    eq_E_str = f"E = {E_const.numerator}/{E_const.denominator} + ({c_E_E9.numerator}/{c_E_E9.denominator}) * E_9"
    if E_const.denominator == 1:
        eq_E_str = f"E = {E_const.numerator} + ({c_E_E9.numerator}/{c_E_E9.denominator}) * E_9"
    print(eq_E_str)
    
    # Print equation for E9
    eq_E9_str = f"E_9 = 1 + ({c_E9_E3_num}/{c_E9_E3_den}) * E_3"
    print(eq_E9_str)
    
    # Print equation for E3
    eq_E3_str = f"E_3 = 1 + ({c_E3_E1_num}/{c_E3_E1_den}) * E_1"
    print(eq_E3_str)
    
    # Print equation for E1
    eq_E1_str = f"E_1 = 2 + ({c_E1_E9_num}/{c_E1_E9_den}) * E_9"
    print(eq_E1_str)
    
    print("\nSolving this system gives:")
    print(f"Minimal expected value of rolls = {E.numerator}/{E.denominator}")

solve_dice_problem()