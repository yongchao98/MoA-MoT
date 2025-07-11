import sympy

def solve_bran_castle_challenge():
    """
    Solves the Bran Castle escape challenge by setting up and solving a system of
    linear equations representing the probabilities of reaching the Treasure Room.
    """
    # Define the symbols for the probabilities of reaching the treasure from each room
    P_ME, P_GH, P_CK, P_SP, P_L, P_TC, P_QR, P_KH = sympy.symbols(
        'P_ME P_GH P_CK P_SP P_L P_TC P_QR P_KH'
    )

    # The probabilities for the absorbing states are known constants.
    # P_TR: Probability from Treasure Room is 1 (Victory)
    # P_VL: Probability from Vampire's Lair is 0 (Game Over)
    P_TR = 1
    P_VL = 0

    # From the castle layout, we derive a system of 8 linear equations.
    # Each equation defines the probability of a room based on its exits.

    # Eq 1: Main Entrance -> 1/2 to Great Hall, 1/2 to Castle Kitchen
    eq1 = sympy.Eq(P_ME, sympy.Rational(1, 2) * P_GH + sympy.Rational(1, 2) * P_CK)

    # Eq 2: Great Hall -> 1/4 to Main Entrance, 1/4 to Secret Passage, 1/4 to Library, 1/4 to Knights' Hall
    eq2 = sympy.Eq(P_GH, sympy.Rational(1, 4) * P_ME + sympy.Rational(1, 4) * P_SP + sympy.Rational(1, 4) * P_L + sympy.Rational(1, 4) * P_KH)

    # Eq 3: Castle Kitchen -> 1/2 to Main Entrance, 1/2 to Knights' Hall
    eq3 = sympy.Eq(P_CK, sympy.Rational(1, 2) * P_ME + sympy.Rational(1, 2) * P_KH)

    # Eq 4: Secret Passage -> 1/4 to Great Hall, 1/4 to Torture Chamber, 1/4 to Treasure Room, 1/4 to Vampire's Lair
    eq4 = sympy.Eq(P_SP, sympy.Rational(1, 4) * P_GH + sympy.Rational(1, 4) * P_TC + sympy.Rational(1, 4) * P_TR + sympy.Rational(1, 4) * P_VL)

    # Eq 5: Library -> 1/3 to Great Hall, 1/3 to Queen's Room, 1/3 to Knights' Hall
    eq5 = sympy.Eq(P_L, sympy.Rational(1, 3) * P_GH + sympy.Rational(1, 3) * P_QR + sympy.Rational(1, 3) * P_KH)

    # Eq 6: Torture Chamber -> 1/2 to Secret Passage, 1/2 to Vampire's Lair
    eq6 = sympy.Eq(P_TC, sympy.Rational(1, 2) * P_SP + sympy.Rational(1, 2) * P_VL)

    # Eq 7: Queen's Room -> 1/3 to Library, 1/3 to Knights' Hall, 1/3 to Treasure Room
    eq7 = sympy.Eq(P_QR, sympy.Rational(1, 3) * P_L + sympy.Rational(1, 3) * P_KH + sympy.Rational(1, 3) * P_TR)

    # Eq 8: Knights' Hall -> 1/4 to Great Hall, 1/4 to Library, 1/4 to Queen's Room, 1/4 to Castle Kitchen
    eq8 = sympy.Eq(P_KH, sympy.Rational(1, 4) * P_GH + sympy.Rational(1, 4) * P_L + sympy.Rational(1, 4) * P_QR + sympy.Rational(1, 4) * P_CK)

    # Solve the system of equations for all unknown probabilities
    solution = sympy.solve(
        [eq1, eq2, eq3, eq4, eq5, eq6, eq7, eq8],
        (P_ME, P_GH, P_CK, P_SP, P_L, P_TC, P_QR, P_KH),
        dict=True
    )[0] # solve returns a list of solutions, we take the first one

    prob_main_entrance = solution[P_ME]
    prob_great_hall = solution[P_GH]
    prob_castle_kitchen = solution[P_CK]

    # Print the final result
    print(f"The exact probability of reaching the Treasure Room from the Main Entrance is {prob_main_entrance}.")
    
    # As requested, show the final numbers in one of the equations to demonstrate the solution.
    print("\nVerification using the equation for the Main Entrance Hall:")
    print(f"P_ME = (1 / 2) * P_GH + (1 / 2) * P_CK")
    print(f"Plugging in the solved values:")
    print(f"{prob_main_entrance} = (1 / 2) * {prob_great_hall} + (1 / 2) * {prob_castle_kitchen}")
    
    # Calculate the right-hand side to show it equals the left-hand side
    rhs = sympy.Rational(1, 2) * prob_great_hall + sympy.Rational(1, 2) * prob_castle_kitchen
    print(f"{prob_main_entrance} = {rhs}")
    
    if prob_main_entrance == rhs:
        print("\nThe equation holds true, confirming the solution.")

solve_bran_castle_challenge()
<<<109/161>>>