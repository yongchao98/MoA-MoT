import sympy

def solve_bran_castle_challenge():
    """
    Calculates the probability of reaching the Treasure Room from the Main Entrance
    by modeling the castle as a system of linear equations and solving it.
    """
    # 1. Define symbolic variables for the probability of reaching the treasure from each room.
    # The subscript corresponds to the room name initials.
    P_ME, P_GH, P_CK, P_SP, P_L, P_KH, P_TC, P_QR = sympy.symbols(
        'P_ME, P_GH, P_CK, P_SP, P_L, P_KH, P_TC, P_QR'
    )
    
    # 2. Set up the system of linear equations based on the castle layout.
    # We use sympy.Rational to maintain exact fractions.
    # P_TR (Treasure Room) = 1, P_VL (Vampire's Lair) = 0.
    
    # Eq 1: Main Entrance -> Great Hall (1/2), Castle Kitchen (1/2)
    eq1 = sympy.Eq(P_ME, sympy.Rational(1, 2) * P_GH + sympy.Rational(1, 2) * P_CK)
    
    # Eq 2: Great Hall -> Main Entrance (1/4), Secret Passage (1/4), Library (1/4), Knights' Hall (1/4)
    eq2 = sympy.Eq(P_GH, sympy.Rational(1, 4) * P_ME + sympy.Rational(1, 4) * P_SP + sympy.Rational(1, 4) * P_L + sympy.Rational(1, 4) * P_KH)
    
    # Eq 3: Secret Passage -> Great Hall (1/4), Torture Chamber (1/4), Treasure Room (1/4), Vampire's Lair (1/4)
    eq3 = sympy.Eq(P_SP, sympy.Rational(1, 4) * P_GH + sympy.Rational(1, 4) * P_TC + sympy.Rational(1, 4) * 1 + sympy.Rational(1, 4) * 0)
    
    # Eq 4: Library -> Great Hall (1/3), Queen's Room (1/3), Knights' Hall (1/3)
    eq4 = sympy.Eq(P_L, sympy.Rational(1, 3) * P_GH + sympy.Rational(1, 3) * P_QR + sympy.Rational(1, 3) * P_KH)
    
    # Eq 5: Torture Chamber -> Secret Passage (1/2), Vampire's Lair (1/2)
    eq5 = sympy.Eq(P_TC, sympy.Rational(1, 2) * P_SP + sympy.Rational(1, 2) * 0)
    
    # Eq 6: Queen's Room -> Library (1/3), Knights' Hall (1/3), Treasure Room (1/3)
    eq6 = sympy.Eq(P_QR, sympy.Rational(1, 3) * P_L + sympy.Rational(1, 3) * P_KH + sympy.Rational(1, 3) * 1)
    
    # Eq 7: Knights' Hall -> Great Hall (1/4), Library (1/4), Queen's Room (1/4), Castle Kitchen (1/4)
    eq7 = sympy.Eq(P_KH, sympy.Rational(1, 4) * P_GH + sympy.Rational(1, 4) * P_L + sympy.Rational(1, 4) * P_QR + sympy.Rational(1, 4) * P_CK)
    
    # Eq 8: Castle Kitchen -> Main Entrance (1/2), Knights' Hall (1/2)
    eq8 = sympy.Eq(P_CK, sympy.Rational(1, 2) * P_ME + sympy.Rational(1, 2) * P_KH)
    
    # 3. Solve the system of equations.
    equations = [eq1, eq2, eq3, eq4, eq5, eq6, eq7, eq8]
    variables = [P_ME, P_GH, P_CK, P_SP, P_L, P_KH, P_TC, P_QR]
    solution = sympy.solve(equations, variables)
    
    # 4. Extract the result for the Main Entrance.
    prob_main_entrance = solution[P_ME]
    
    # 5. Print the final answer as an irreducible fraction.
    numerator = prob_main_entrance.p
    denominator = prob_main_entrance.q
    
    print("The exact probability of reaching the Treasure Room from the Main Entrance is given by the equation:")
    print(f"P(Success) = {numerator} / {denominator}")

if __name__ == '__main__':
    solve_bran_castle_challenge()