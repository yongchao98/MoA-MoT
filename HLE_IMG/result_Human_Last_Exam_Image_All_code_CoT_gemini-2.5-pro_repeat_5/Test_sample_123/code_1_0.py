import sympy

def solve_castle_escape():
    """
    Calculates the exact probability of reaching the Treasure Room from the Main Entrance
    by solving a system of linear equations representing the castle as a Markov chain.
    """
    # Define the symbols for the probabilities of reaching the Treasure Room from each room
    P_ME, P_GH, P_CK, P_SP, P_L, P_KH, P_TC, P_QR = sympy.symbols(
        'P_ME P_GH P_CK P_SP P_L P_KH P_TC P_QR'
    )

    # Define the system of linear equations based on the transition probabilities.
    # We use sympy.Rational for exact fractional arithmetic.
    # The probabilities for the absorbing states are P_TreasureRoom = 1 and P_VampireLair = 0.

    # Eq 1: Main Entrance -> Great Hall (1/2), Castle Kitchen (1/2)
    eq1 = sympy.Eq(P_ME, sympy.Rational(1, 2) * P_GH + sympy.Rational(1, 2) * P_CK)

    # Eq 2: Great Hall -> Main Entrance (1/4), Secret Passage (1/4), Library (1/4), Knights' Hall (1/4)
    eq2 = sympy.Eq(P_GH, sympy.Rational(1, 4) * P_ME + sympy.Rational(1, 4) * P_SP + sympy.Rational(1, 4) * P_L + sympy.Rational(1, 4) * P_KH)

    # Eq 3: Castle Kitchen -> Main Entrance (1/2), Knights' Hall (1/2)
    eq3 = sympy.Eq(P_CK, sympy.Rational(1, 2) * P_ME + sympy.Rational(1, 2) * P_KH)

    # Eq 4: Secret Passage -> Great Hall (1/4), Torture Chamber (1/4), Treasure Room (1/4), Vampire's Lair (1/4)
    eq4 = sympy.Eq(P_SP, sympy.Rational(1, 4) * P_GH + sympy.Rational(1, 4) * P_TC + sympy.Rational(1, 4) * 1)

    # Eq 5: Library -> Great Hall (1/3), Queen's Room (1/3), Knights' Hall (1/3)
    eq5 = sympy.Eq(P_L, sympy.Rational(1, 3) * P_GH + sympy.Rational(1, 3) * P_QR + sympy.Rational(1, 3) * P_KH)

    # Eq 6: Knights' Hall -> Great Hall (1/4), Library (1/4), Queen's Room (1/4), Castle Kitchen (1/4)
    eq6 = sympy.Eq(P_KH, sympy.Rational(1, 4) * P_GH + sympy.Rational(1, 4) * P_L + sympy.Rational(1, 4) * P_QR + sympy.Rational(1, 4) * P_CK)

    # Eq 7: Torture Chamber -> Secret Passage (1/2), Vampire's Lair (1/2)
    eq7 = sympy.Eq(P_TC, sympy.Rational(1, 2) * P_SP)

    # Eq 8: Queen's Room -> Library (1/3), Knights' Hall (1/3), Treasure Room (1/3)
    eq8 = sympy.Eq(P_QR, sympy.Rational(1, 3) * P_L + sympy.Rational(1, 3) * P_KH + sympy.Rational(1, 3) * 1)

    # Solve the system of 8 linear equations for the 8 unknown probabilities
    solution = sympy.solve(
        (eq1, eq2, eq3, eq4, eq5, eq6, eq7, eq8),
        (P_ME, P_GH, P_CK, P_SP, P_L, P_KH, P_TC, P_QR)
    )

    # The probability we are looking for is for the Main Entrance
    prob_main_entrance = solution[P_ME]
    
    # Print the final result as an irreducible fraction
    numerator = prob_main_entrance.p
    denominator = prob_main_entrance.q
    
    print(f"The system of equations is solved to find the probability of reaching the Treasure Room.")
    print(f"The probability of success starting from the Main Entrance is P_ME.")
    print(f"The final equation for the probability is: P_ME = {numerator} / {denominator}")

if __name__ == "__main__":
    solve_castle_escape()