import sympy

def solve_bran_castle_challenge():
    """
    Calculates the exact probability of reaching the Treasure Room from the Main Entrance
    by setting up and solving a system of linear equations.
    """
    # 1. Define the variables for the probabilities of reaching the Treasure Room from each room.
    # The probabilities for the absorbing states are known: P_TR = 1, P_VL = 0.
    P_ME, P_GH, P_CK, P_KH, P_L, P_QR, P_SP, P_TC = sympy.symbols(
        'P_ME P_GH P_CK P_KH P_L P_QR P_SP P_TC'
    )

    # 2. Set up the system of 8 linear equations based on the castle layout and transition probabilities.
    # Each equation represents the probability of success from one room.
    # We use sympy.Rational to ensure exact fractional arithmetic.
    
    # Equation for Main Entrance
    eq1 = sympy.Eq(P_ME, sympy.Rational(1, 2) * P_GH + sympy.Rational(1, 2) * P_CK)
    
    # Equation for Great Hall
    eq2 = sympy.Eq(P_GH, sympy.Rational(1, 4) * P_ME + sympy.Rational(1, 4) * P_SP + sympy.Rational(1, 4) * P_L + sympy.Rational(1, 4) * P_KH)
    
    # Equation for Secret Passage (P_TR=1, P_VL=0)
    eq3 = sympy.Eq(P_SP, sympy.Rational(1, 4) * P_GH + sympy.Rational(1, 4) * P_TC + sympy.Rational(1, 4) * 1 + sympy.Rational(1, 4) * 0)
    
    # Equation for Library
    eq4 = sympy.Eq(P_L, sympy.Rational(1, 3) * P_GH + sympy.Rational(1, 3) * P_QR + sympy.Rational(1, 3) * P_KH)
    
    # Equation for Torture Chamber (P_VL=0)
    eq5 = sympy.Eq(P_TC, sympy.Rational(1, 2) * P_SP + sympy.Rational(1, 2) * 0)
    
    # Equation for Queen's Room (P_TR=1)
    eq6 = sympy.Eq(P_QR, sympy.Rational(1, 3) * P_L + sympy.Rational(1, 3) * P_KH + sympy.Rational(1, 3) * 1)
    
    # Equation for Knights' Hall
    eq7 = sympy.Eq(P_KH, sympy.Rational(1, 4) * P_GH + sympy.Rational(1, 4) * P_L + sympy.Rational(1, 4) * P_QR + sympy.Rational(1, 4) * P_CK)
    
    # Equation for Castle Kitchen
    eq8 = sympy.Eq(P_CK, sympy.Rational(1, 2) * P_ME + sympy.Rational(1, 2) * P_KH)

    # 3. Solve the system of equations for all unknown probabilities.
    # The solution will be a dictionary mapping each variable to its value.
    solution = sympy.solve(
        [eq1, eq2, eq3, eq4, eq5, eq6, eq7, eq8],
        (P_ME, P_GH, P_CK, P_KH, P_L, P_QR, P_SP, P_TC)
    )

    # 4. Extract the probability for the Main Entrance and present the result.
    p_main_entrance = solution[P_ME]
    
    numerator, denominator = p_main_entrance.as_numer_denom()

    print("The system of linear equations representing the problem has been solved.")
    print(f"The probability of eventually reaching the Treasure Room from the Main Entrance is:")
    print(f"P(Main Entrance) = {numerator} / {denominator}")

if __name__ == '__main__':
    solve_bran_castle_challenge()