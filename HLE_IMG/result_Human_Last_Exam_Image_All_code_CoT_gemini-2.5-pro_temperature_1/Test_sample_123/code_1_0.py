import sympy

def solve_bran_castle_escape():
    """
    Calculates the exact probability of reaching the Treasure Room
    from the Main Entrance by setting up and solving a system of linear equations.
    """
    # Define the symbols for the probability of reaching the treasure from each room
    P_ME, P_GH, P_CK, P_SP, P_L, P_KH, P_TC, P_QR = sympy.symbols(
        'P_ME P_GH P_CK P_SP P_L P_KH P_TC P_QR'
    )

    # Set up the system of linear equations based on the castle layout
    # The probability of reaching the goal from a room is the sum of
    # (prob_to_move_to_next_room * prob_of_reaching_goal_from_next_room)
    # P_TreasureRoom = 1, P_VampireLair = 0

    # Equation for Main Entrance (2 exits)
    # Exits: Great Hall, Castle Kitchen
    eq_ME = sympy.Eq(P_ME, sympy.Rational(1, 2) * P_GH + sympy.Rational(1, 2) * P_CK)
    
    # Equation for Great Hall (4 exits)
    # Exits: Main Entrance, Secret Passage, Library, Knights' Hall
    eq_GH = sympy.Eq(P_GH, sympy.Rational(1, 4) * P_ME + sympy.Rational(1, 4) * P_SP + sympy.Rational(1, 4) * P_L + sympy.Rational(1, 4) * P_KH)

    # Equation for Castle Kitchen (2 exits)
    # Exits: Main Entrance, Knights' Hall
    eq_CK = sympy.Eq(P_CK, sympy.Rational(1, 2) * P_ME + sympy.Rational(1, 2) * P_KH)

    # Equation for Secret Passage (4 exits)
    # Exits: Great Hall, Torture Chamber, Treasure Room (Prob=1), Vampire's Lair (Prob=0)
    eq_SP = sympy.Eq(P_SP, sympy.Rational(1, 4) * P_GH + sympy.Rational(1, 4) * P_TC + sympy.Rational(1, 4) * 1)

    # Equation for Library (3 exits)
    # Exits: Great Hall, Queen's Room, Knights' Hall
    eq_L = sympy.Eq(P_L, sympy.Rational(1, 3) * P_GH + sympy.Rational(1, 3) * P_QR + sympy.Rational(1, 3) * P_KH)
    
    # Equation for Knights' Hall (4 exits)
    # Exits: Great Hall, Library, Queen's Room, Castle Kitchen
    eq_KH = sympy.Eq(P_KH, sympy.Rational(1, 4) * P_GH + sympy.Rational(1, 4) * P_L + sympy.Rational(1, 4) * P_QR + sympy.Rational(1, 4) * P_CK)

    # Equation for Torture Chamber (2 exits)
    # Exits: Secret Passage, Vampire's Lair (Prob=0)
    eq_TC = sympy.Eq(P_TC, sympy.Rational(1, 2) * P_SP)
    
    # Equation for Queen's Room (3 exits)
    # Exits: Library, Knights' Hall, Treasure Room (Prob=1)
    eq_QR = sympy.Eq(P_QR, sympy.Rational(1, 3) * P_L + sympy.Rational(1, 3) * P_KH + sympy.Rational(1, 3) * 1)

    # Print the equations for clarity
    print("System of Linear Equations:")
    print(f"1. P_ME = {eq_ME.rhs}")
    print(f"2. P_GH = {eq_GH.rhs}")
    print(f"3. P_CK = {eq_CK.rhs}")
    print(f"4. P_SP = {eq_SP.rhs}")
    print(f"5. P_L  = {eq_L.rhs}")
    print(f"6. P_KH = {eq_KH.rhs}")
    print(f"7. P_TC = {eq_TC.rhs}")
    print(f"8. P_QR = {eq_QR.rhs}")
    
    # Solve the system of equations
    solution = sympy.solve(
        [eq_ME, eq_GH, eq_CK, eq_SP, eq_L, eq_KH, eq_TC, eq_QR],
        [P_ME, P_GH, P_CK, P_SP, P_L, P_KH, P_TC, P_QR]
    )

    # Extract and print the final probability
    final_prob = solution[P_ME]
    print("\nSolving this system gives the following probability for the Main Entrance:")
    print(f"P_ME = {final_prob.p}/{final_prob.q}")

solve_bran_castle_escape()
<<<109/161>>>