import sympy

def solve_bran_castle_escape():
    """
    Calculates the exact probability of escaping with the treasure from Bran Castle.
    """
    # Define symbols for the unknown probabilities of reaching the Treasure Room from each room.
    # The name corresponds to the room's initials.
    P_ME, P_GH, P_CK, P_SP, P_L, P_KH, P_QR, P_TC = sympy.symbols(
        'P_ME P_GH P_CK P_SP P_L P_KH P_QR P_TC'
    )

    # Define the probabilities for the absorbing states.
    # P_TR: Probability of success from the Treasure Room is 1.
    # P_VL: Probability of success from the Vampire's Lair is 0.
    P_TR = 1
    P_VL = 0

    # The problem states that from any room, each exit has an equal probability of being chosen.
    # We set up a system of linear equations based on the connections described.
    # The probability of success from a room is the sum of (prob of taking an exit * prob of success from the destination).
    
    # Equation for Main Entrance (2 exits)
    eq1 = sympy.Eq(P_ME, sympy.Rational(1, 2) * P_GH + sympy.Rational(1, 2) * P_CK)
    
    # Equation for Great Hall (4 exits)
    eq2 = sympy.Eq(P_GH, sympy.Rational(1, 4) * P_ME + sympy.Rational(1, 4) * P_SP + sympy.Rational(1, 4) * P_L + sympy.Rational(1, 4) * P_KH)
    
    # Equation for Castle Kitchen (2 exits)
    eq3 = sympy.Eq(P_CK, sympy.Rational(1, 2) * P_ME + sympy.Rational(1, 2) * P_KH)
    
    # Equation for Secret Passage (4 exits, including absorbing states)
    eq4 = sympy.Eq(P_SP, sympy.Rational(1, 4) * P_GH + sympy.Rational(1, 4) * P_TC + sympy.Rational(1, 4) * P_TR + sympy.Rational(1, 4) * P_VL)
    
    # Equation for Library (3 exits)
    eq5 = sympy.Eq(P_L, sympy.Rational(1, 3) * P_GH + sympy.Rational(1, 3) * P_QR + sympy.Rational(1, 3) * P_KH)
    
    # Equation for Knights' Hall (4 exits)
    eq6 = sympy.Eq(P_KH, sympy.Rational(1, 4) * P_GH + sympy.Rational(1, 4) * P_L + sympy.Rational(1, 4) * P_QR + sympy.Rational(1, 4) * P_CK)
    
    # Equation for Queen's Room (3 exits, including an absorbing state)
    eq7 = sympy.Eq(P_QR, sympy.Rational(1, 3) * P_L + sympy.Rational(1, 3) * P_KH + sympy.Rational(1, 3) * P_TR)
    
    # Equation for Torture Chamber (2 exits, including an absorbing state)
    eq8 = sympy.Eq(P_TC, sympy.Rational(1, 2) * P_SP + sympy.Rational(1, 2) * P_VL)

    # Use sympy's solve function to find the values of all unknown probabilities.
    solution = sympy.solve(
        [eq1, eq2, eq3, eq4, eq5, eq6, eq7, eq8],
        (P_ME, P_GH, P_CK, P_SP, P_L, P_KH, P_QR, P_TC),
        dict=True
    )[0] # solve returns a list of solutions, we take the first one

    # The result we want is the probability of success starting from the Main Entrance (P_ME).
    prob_main_entrance = solution[P_ME]
    
    # Extract the numerator and denominator for printing.
    numerator = prob_main_entrance.p
    denominator = prob_main_entrance.q

    print("The probability of reaching the Treasure Room from the Main Entrance is P_ME.")
    print("Solving the system of linear equations for P_ME gives the following result:")
    print(f"P_ME = {numerator} / {denominator}")

solve_bran_castle_escape()
<<<109/161>>>