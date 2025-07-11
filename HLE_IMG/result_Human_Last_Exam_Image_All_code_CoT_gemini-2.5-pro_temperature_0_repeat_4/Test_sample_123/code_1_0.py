import sympy

def solve_castle_escape():
    """
    Sets up and solves the system of linear equations for the Bran Castle escape problem.
    """
    # Define the symbols for the unknown probabilities for each room
    p_ME, p_GH, p_CK, p_SP, p_L, p_KH, p_TC, p_QR = sympy.symbols(
        'p_ME p_GH p_CK p_SP p_L p_KH p_TC p_QR'
    )

    # Use sympy.Rational for exact fractions
    R = sympy.Rational

    # Define the system of linear equations based on the transition probabilities.
    # We rearrange each equation to the form f(p) = 0.
    # e.g., p_ME = (1/2)p_GH + (1/2)p_CK  becomes  p_ME - (1/2)p_GH - (1/2)p_CK = 0
    
    # Equation for Main Entrance
    eq1 = sympy.Eq(p_ME - R(1, 2) * p_GH - R(1, 2) * p_CK, 0)
    
    # Equation for Great Hall
    eq2 = sympy.Eq(p_GH - R(1, 4) * p_ME - R(1, 4) * p_SP - R(1, 4) * p_L - R(1, 4) * p_KH, 0)
    
    # Equation for Castle Kitchen
    eq3 = sympy.Eq(p_CK - R(1, 2) * p_ME - R(1, 2) * p_KH, 0)
    
    # Equation for Secret Passage (p_Treasure=1, p_Vampire=0)
    eq4 = sympy.Eq(p_SP - R(1, 4) * p_GH - R(1, 4) * p_TC, R(1, 4) * 1 + R(1, 4) * 0)
    
    # Equation for Library
    eq5 = sympy.Eq(p_L - R(1, 3) * p_GH - R(1, 3) * p_QR - R(1, 3) * p_KH, 0)
    
    # Equation for Knights' Hall
    eq6 = sympy.Eq(p_KH - R(1, 4) * p_GH - R(1, 4) * p_L - R(1, 4) * p_QR - R(1, 4) * p_CK, 0)
    
    # Equation for Torture Chamber (p_Vampire=0)
    eq7 = sympy.Eq(p_TC - R(1, 2) * p_SP, R(1, 2) * 0)
    
    # Equation for Queen's Room (p_Treasure=1)
    eq8 = sympy.Eq(p_QR - R(1, 3) * p_L - R(1, 3) * p_KH, R(1, 3) * 1)

    equations = [eq1, eq2, eq3, eq4, eq5, eq6, eq7, eq8]
    
    print("System of Linear Equations to be Solved:")
    print("p_Treasure = 1")
    print("p_Vampire = 0")
    # The instruction "output each number in the final equation" is interpreted as showing the system being solved.
    print(f"1. p_ME = {R(1,2)}*p_GH + {R(1,2)}*p_CK")
    print(f"2. p_GH = {R(1,4)}*p_ME + {R(1,4)}*p_SP + {R(1,4)}*p_L + {R(1,4)}*p_KH")
    print(f"3. p_CK = {R(1,2)}*p_ME + {R(1,2)}*p_KH")
    print(f"4. p_SP = {R(1,4)}*p_GH + {R(1,4)}*p_TC + {R(1,4)}")
    print(f"5. p_L = {R(1,3)}*p_GH + {R(1,3)}*p_QR + {R(1,3)}*p_KH")
    print(f"6. p_KH = {R(1,4)}*p_GH + {R(1,4)}*p_L + {R(1,4)}*p_QR + {R(1,4)}*p_CK")
    print(f"7. p_TC = {R(1,2)}*p_SP")
    print(f"8. p_QR = {R(1,3)}*p_L + {R(1,3)}*p_KH + {R(1,3)}")
    print("-" * 20)

    # Solve the system of equations
    solution = sympy.solve(equations, (p_ME, p_GH, p_CK, p_SP, p_L, p_KH, p_TC, p_QR))

    # Extract the probability for the Main Entrance
    prob_main_entrance = solution[p_ME]
    
    # The final equation is the solved value for p_ME
    print(f"The final calculated probability for the Main Entrance is:")
    print(f"p_ME = {prob_main_entrance.p} / {prob_main_entrance.q}")


solve_castle_escape()