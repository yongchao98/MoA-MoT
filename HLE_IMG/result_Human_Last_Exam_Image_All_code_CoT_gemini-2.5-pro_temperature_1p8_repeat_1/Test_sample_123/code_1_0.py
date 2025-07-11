import sympy
from fractions import Fraction

def solve_castle_escape():
    """
    Solves the Bran Castle escape problem by setting up and solving a system of linear equations.
    """
    # 1. Define the variables for the probability of success from each room.
    P_ME, P_GH, P_CK, P_SP, P_L, P_TC, P_QR, P_KH = sympy.symbols(
        'P_ME P_GH P_CK P_SP P_L P_TC P_QR P_KH'
    )

    # The probabilities from the absorbing states are known constants.
    P_TR = 1  # Treasure Room (Success)
    P_VL = 0  # Vampire's Lair (Failure)

    # 2. Set up the system of linear equations based on the graph transitions.
    # Note: Using sympy.Rational ensures calculations are exact fractions.
    eq1 = sympy.Eq(P_ME, sympy.Rational(1, 2) * P_GH + sympy.Rational(1, 2) * P_CK)
    eq2 = sympy.Eq(P_GH, sympy.Rational(1, 4) * P_ME + sympy.Rational(1, 4) * P_SP + sympy.Rational(1, 4) * P_L + sympy.Rational(1, 4) * P_KH)
    eq3 = sympy.Eq(P_CK, sympy.Rational(1, 2) * P_ME + sympy.Rational(1, 2) * P_KH)
    eq4 = sympy.Eq(P_SP, sympy.Rational(1, 4) * P_GH + sympy.Rational(1, 4) * P_TC + sympy.Rational(1, 4) * P_TR)
    eq5 = sympy.Eq(P_L, sympy.Rational(1, 3) * P_GH + sympy.Rational(1, 3) * P_QR + sympy.Rational(1, 3) * P_KH)
    eq6 = sympy.Eq(P_TC, sympy.Rational(1, 2) * P_SP + sympy.Rational(1, 2) * P_VL)
    eq7 = sympy.Eq(P_QR, sympy.Rational(1, 3) * P_L + sympy.Rational(1, 3) * P_KH + sympy.Rational(1, 3) * P_TR)
    eq8 = sympy.Eq(P_KH, sympy.Rational(1, 4) * P_GH + sympy.Rational(1, 4) * P_L + sympy.Rational(1, 4) * P_QR + sympy.Rational(1, 4) * P_CK)

    # 3. Solve the system of equations.
    equations = [eq1, eq2, eq3, eq4, eq5, eq6, eq7, eq8]
    variables = [P_ME, P_GH, P_CK, P_SP, P_L, P_TC, P_QR, P_KH]
    solution = sympy.solve(equations, variables)

    # 4. Extract the probabilities needed for the final calculation display.
    prob_me = solution[P_ME]
    prob_gh = solution[P_GH]
    prob_ck = solution[P_CK]
    
    # Use Fraction for clean numerator/denominator representation
    frac_me = Fraction(prob_me).limit_denominator()
    frac_gh = Fraction(prob_gh).limit_denominator()
    frac_ck = Fraction(prob_ck).limit_denominator()

    # 5. Print the results and the final equation steps.
    print("The final probability is calculated from the initial room, the Main Entrance.")
    print("The equation for the Main Entrance is:")
    print(f"P_ME = 1/2 * P_GH + 1/2 * P_CK")
    print("\nAfter solving the full system of equations, we find the probabilities for the connecting rooms:")
    print(f"P_GH = {frac_gh.numerator}/{frac_gh.denominator}")
    print(f"P_CK = {frac_ck.numerator}/{frac_ck.denominator}")
    print("\nSubstituting these into the equation for P_ME:")
    print(f"P_ME = (1/2) * ({frac_gh.numerator}/{frac_gh.denominator}) + (1/2) * ({frac_ck.numerator}/{frac_ck.denominator})")
    
    # Calculate the combined fraction to show the final steps
    combined_num = frac_gh.numerator * frac_ck.denominator + frac_ck.numerator * frac_gh.denominator
    combined_den = 2 * frac_gh.denominator * frac_ck.denominator
    final_frac = Fraction(combined_num, combined_den).limit_denominator()

    print(f"P_ME = {frac_gh.numerator}/{2 * frac_gh.denominator} + {frac_ck.numerator}/{2 * frac_ck.denominator}")
    print(f"P_ME = {final_frac.numerator}/{final_frac.denominator}")
    print("\nThe final probability of reaching the Treasure Room is:")
    print(f"{frac_me.numerator}/{frac_me.denominator}")

solve_castle_escape()