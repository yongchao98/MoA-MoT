import sympy

def solve_bran_castle_challenge():
    """
    This function models the Bran Castle escape challenge as a system of linear
    equations and solves it to find the exact probability of success.
    """
    # Define symbolic variables for the probability of success from each room.
    # The variable name corresponds to the first letter(s) of the room name.
    # ME: Main Entrance, GH: Great Hall, CK: Castle Kitchen, SP: Secret Passage,
    # L: Library, KH: Knights' Hall, TC: Torture Chamber, QR: Queen's Room.
    P_ME, P_GH, P_CK, P_SP, P_L, P_KH, P_TC, P_QR = sympy.symbols(
        'P_ME, P_GH, P_CK, P_SP, P_L, P_KH, P_TC, P_QR'
    )

    # Probabilities for the absorbing states (the end points).
    P_TreasureRoom = 1
    P_VampireLair = 0

    # The system of linear equations is derived from the castle layout.
    # For each room, P(Room) = sum[ P(NextRoom) * (1 / num_exits) ]
    equations = [
        # Main Entrance -> Great Hall (1/2), Castle Kitchen (1/2)
        sympy.Eq(P_ME, sympy.Rational(1, 2) * P_GH + sympy.Rational(1, 2) * P_CK),

        # Great Hall -> Main Entrance (1/4), Secret Passage (1/4), Library (1/4), Knights' Hall (1/4)
        sympy.Eq(P_GH, sympy.Rational(1, 4) * P_ME + sympy.Rational(1, 4) * P_SP + sympy.Rational(1, 4) * P_L + sympy.Rational(1, 4) * P_KH),

        # Castle Kitchen -> Main Entrance (1/2), Knights' Hall (1/2)
        sympy.Eq(P_CK, sympy.Rational(1, 2) * P_ME + sympy.Rational(1, 2) * P_KH),

        # Secret Passage -> Great Hall (1/4), Torture Chamber (1/4), Treasure Room (1/4), Vampire's Lair (1/4)
        sympy.Eq(P_SP, sympy.Rational(1, 4) * P_GH + sympy.Rational(1, 4) * P_TC + sympy.Rational(1, 4) * P_TreasureRoom + sympy.Rational(1, 4) * P_VampireLair),

        # Library -> Great Hall (1/3), Queen's Room (1/3), Knights' Hall (1/3)
        sympy.Eq(P_L, sympy.Rational(1, 3) * P_GH + sympy.Rational(1, 3) * P_QR + sympy.Rational(1, 3) * P_KH),

        # Knights' Hall -> Great Hall (1/4), Library (1/4), Queen's Room (1/4), Castle Kitchen (1/4)
        sympy.Eq(P_KH, sympy.Rational(1, 4) * P_GH + sympy.Rational(1, 4) * P_L + sympy.Rational(1, 4) * P_QR + sympy.Rational(1, 4) * P_CK),

        # Torture Chamber -> Secret Passage (1/2), Vampire's Lair (1/2)
        sympy.Eq(P_TC, sympy.Rational(1, 2) * P_SP + sympy.Rational(1, 2) * P_VampireLair),

        # Queen's Room -> Library (1/3), Knights' Hall (1/3), Treasure Room (1/3)
        sympy.Eq(P_QR, sympy.Rational(1, 3) * P_L + sympy.Rational(1, 3) * P_KH + sympy.Rational(1, 3) * P_TreasureRoom)
    ]

    # Solve the system of equations for all unknown probabilities.
    solution = sympy.solve(equations, (P_ME, P_GH, P_CK, P_SP, P_L, P_KH, P_TC, P_QR))

    # The required probability is for the Main Entrance.
    prob_main_entrance = solution[P_ME]
    
    # Extract numerator and denominator for the final equation.
    numerator = prob_main_entrance.p
    denominator = prob_main_entrance.q

    print(f"The probability of reaching the Treasure Room from the Main Entrance is {prob_main_entrance}.")
    print("\nThis result is derived from solving a system of 8 linear equations.")
    print("The final probability P_ME can be represented by the equation:")
    print(f"{denominator} * P_ME = {numerator}")

solve_bran_castle_challenge()
<<<109/161>>>