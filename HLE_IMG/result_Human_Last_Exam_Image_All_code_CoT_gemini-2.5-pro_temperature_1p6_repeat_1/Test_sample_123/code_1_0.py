import sympy

def solve_bran_castle_escape():
    """
    Solves the Bran Castle escape probability problem using a system of linear equations.
    """
    # Define the symbols for the probabilities of reaching the Treasure Room from each room.
    # P_ME: Main Entrance, P_GH: Great Hall, P_CK: Castle Kitchen, P_SP: Secret Passage,
    # P_L: Library, P_KH: Knights' Hall, P_TC: Torture Chamber, P_QR: Queen's Room.
    P_ME, P_GH, P_CK, P_SP, P_L, P_KH, P_TC, P_QR = sympy.symbols(
        'P_ME P_GH P_CK P_SP P_L P_KH P_TC P_QR'
    )

    # Probabilities for the absorbing states
    # P_TR: Treasure Room (Victory)
    P_TR = 1
    # P_VL: Vampire's Lair (Game Over)
    P_VL = 0

    # System of linear equations based on the problem description.
    # Each exit from a room has equal probability.

    # Equation 1: Main Entrance (2 exits: Great Hall, Castle Kitchen)
    eq1 = sympy.Eq(P_ME, sympy.S(1)/2 * P_GH + sympy.S(1)/2 * P_CK)

    # Equation 2: Great Hall (4 exits: Main Entrance, Secret Passage, Library, Knights' Hall)
    eq2 = sympy.Eq(P_GH, sympy.S(1)/4 * P_ME + sympy.S(1)/4 * P_SP + sympy.S(1)/4 * P_L + sympy.S(1)/4 * P_KH)

    # Equation 3: Castle Kitchen (2 exits: Main Entrance, Knights' Hall)
    eq3 = sympy.Eq(P_CK, sympy.S(1)/2 * P_ME + sympy.S(1)/2 * P_KH)

    # Equation 4: Secret Passage (4 exits: Great Hall, Torture Chamber, Treasure Room, Vampire's Lair)
    eq4 = sympy.Eq(P_SP, sympy.S(1)/4 * P_GH + sympy.S(1)/4 * P_TC + sympy.S(1)/4 * P_TR + sympy.S(1)/4 * P_VL)

    # Equation 5: Library (3 exits: Great Hall, Queen's Room, Knights' Hall)
    eq5 = sympy.Eq(P_L, sympy.S(1)/3 * P_GH + sympy.S(1)/3 * P_QR + sympy.S(1)/3 * P_KH)
    
    # Equation 6: Knights' Hall (4 exits: Great Hall, Library, Queen's Room, Castle Kitchen)
    eq6 = sympy.Eq(P_KH, sympy.S(1)/4 * P_GH + sympy.S(1)/4 * P_L + sympy.S(1)/4 * P_QR + sympy.S(1)/4 * P_CK)

    # Equation 7: Torture Chamber (2 exits: Secret Passage, Vampire's Lair)
    eq7 = sympy.Eq(P_TC, sympy.S(1)/2 * P_SP + sympy.S(1)/2 * P_VL)

    # Equation 8: Queen Marie's Room (3 exits: Library, Knights' Hall, Treasure Room)
    eq8 = sympy.Eq(P_QR, sympy.S(1)/3 * P_L + sympy.S(1)/3 * P_KH + sympy.S(1)/3 * P_TR)
    
    print("This problem can be modeled as a system of linear equations, where each variable represents the probability of reaching the Treasure Room from a specific location in the castle.")
    print("Let P_X be the probability of success starting from room X.")
    print("\nThe absorbing states probabilities are:")
    print("P_TreasureRoom = 1")
    print("P_VampireLair = 0\n")
    print("The system of equations for the other rooms is as follows:")
    print(f"1. {eq1}")
    print(f"2. {eq2}")
    print(f"3. {eq3}")
    print(f"4. {eq4.subs([(P_TR, 1), (P_VL, 0)])}")
    print(f"5. {eq5}")
    print(f"6. {eq6}")
    print(f"7. {eq7.subs(P_VL, 0)}")
    print(f"8. {eq8.subs(P_TR, 1)}")
    
    # Solve the system of equations
    solution = sympy.solve([eq1, eq2, eq3, eq4, eq5, eq6, eq7, eq8],
                           (P_ME, P_GH, P_CK, P_SP, P_L, P_KH, P_TC, P_QR))

    # The probability of success from the Main Entrance is the value of P_ME
    prob_me = solution[P_ME]
    
    # Extract numerator and denominator
    numerator = prob_me.p
    denominator = prob_me.q

    print("\nSolving this system for P_ME gives the probability of reaching the Treasure Room from the Main Entrance.")
    print("The exact probability is given as an irreducible fraction.")
    print("\nThe final equation for the probability from the Main Entrance is:")
    print(f"P_ME = {numerator} / {denominator}")

# Execute the function
solve_bran_castle_escape()