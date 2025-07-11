import sympy

def solve_bran_castle_challenge():
    """
    Calculates the exact probability of reaching the Treasure Room from the Main Entrance
    by modeling the castle as a system of linear equations and solving it using sympy.
    """
    # Define symbols for the unknown probabilities for each non-absorbing room
    p_ME, p_GH, p_SP, p_L, p_TC, p_QR, p_KH, p_CK = sympy.symbols(
        'p_ME p_GH p_SP p_L p_TC p_QR p_KH p_CK'
    )

    # Use sympy.Rational for exact fractions to avoid floating point errors.
    # Alias sympy.Rational as F for convenience.
    F = sympy.Rational

    # Set up the system of 8 linear equations based on the rules.
    # The probability of winning from a room is the weighted average of the probabilities
    # of winning from the rooms it connects to.
    # The weights are the transition probabilities.
    # P_TreasureRoom = 1, P_VampiresLair = 0

    # 1. Main Entrance (connects to Great Hall, Castle Kitchen)
    eq1 = sympy.Eq(p_ME, F(1, 2) * p_GH + F(1, 2) * p_CK)

    # 2. Great Hall (connects to Main Entrance, Secret Passage, Library, Knights' Hall)
    eq2 = sympy.Eq(p_GH, F(1, 4) * p_ME + F(1, 4) * p_SP + F(1, 4) * p_L + F(1, 4) * p_KH)

    # 3. Secret Passage (connects to Great Hall, Torture Chamber, Treasure Room, Vampire's Lair)
    eq3 = sympy.Eq(p_SP, F(1, 4) * p_GH + F(1, 4) * p_TC + F(1, 4) * 1 + F(1, 4) * 0)

    # 4. Library (connects to Great Hall, Queen's Room, Knights' Hall)
    eq4 = sympy.Eq(p_L, F(1, 3) * p_GH + F(1, 3) * p_QR + F(1, 3) * p_KH)

    # 5. Torture Chamber (connects to Secret Passage, Vampire's Lair)
    eq5 = sympy.Eq(p_TC, F(1, 2) * p_SP + F(1, 2) * 0)

    # 6. Queen's Room (connects to Library, Knights' Hall, Treasure Room)
    eq6 = sympy.Eq(p_QR, F(1, 3) * p_L + F(1, 3) * p_KH + F(1, 3) * 1)

    # 7. Knights' Hall (connects to Great Hall, Library, Queen's Room, Castle Kitchen)
    eq7 = sympy.Eq(p_KH, F(1, 4) * p_GH + F(1, 4) * p_L + F(1, 4) * p_QR + F(1, 4) * p_CK)

    # 8. Castle Kitchen (connects to Main Entrance, Knights' Hall)
    eq8 = sympy.Eq(p_CK, F(1, 2) * p_ME + F(1, 2) * p_KH)

    # Solve the system of 8 equations for the 8 unknown probabilities
    equations = [eq1, eq2, eq3, eq4, eq5, eq6, eq7, eq8]
    variables = [p_ME, p_GH, p_SP, p_L, p_TC, p_QR, p_KH, p_CK]
    solution = sympy.solve(equations, variables)

    # The result we want is the probability from the Main Entrance
    result_fraction = solution[p_ME]

    # Print the result in a descriptive format, including the numerator and denominator
    numerator = result_fraction.p
    denominator = result_fraction.q
    
    print(f"The exact probability of reaching the Treasure Room from the Main Entrance is {result_fraction}.")
    print(f"The equation is: Probability = {numerator} / {denominator}")

if __name__ == '__main__':
    solve_bran_castle_challenge()