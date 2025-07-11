import sympy

def solve_castle_escape():
    """
    Sets up and solves the system of linear equations for the Bran Castle escape problem
    to find the exact probability of reaching the Treasure Room from the Main Entrance.
    """
    # Define the variables (probabilities for each room)
    # p0: Main Entrance, p1: Great Hall, p2: Castle Kitchen, p3: Secret Passage,
    # p4: Library, p5: Knights' Hall, p6: Torture Chamber, p7: Queen's Room
    p0, p1, p2, p3, p4, p5, p6, p7 = sympy.symbols('p0 p1 p2 p3 p4 p5 p6 p7')

    # Probabilities for absorbing states
    p_treasure = 1
    p_vampire = 0

    # Define the system of linear equations based on the transition probabilities
    # sympy.Eq(LHS, RHS) creates an equation object. Fractions are kept exact.
    eq1 = sympy.Eq(p0, p1/2 + p2/2)
    eq2 = sympy.Eq(p1, p0/4 + p3/4 + p4/4 + p5/4)
    eq3 = sympy.Eq(p2, p0/2 + p5/2)
    eq4 = sympy.Eq(p3, p1/4 + p6/4 + p_treasure/4 + p_vampire/4)
    eq5 = sympy.Eq(p4, p1/3 + p7/3 + p5/3)
    eq6 = sympy.Eq(p5, p1/4 + p2/4 + p4/4 + p7/4)
    eq7 = sympy.Eq(p6, p3/2 + p_vampire/2)
    eq8 = sympy.Eq(p7, p4/3 + p5/3 + p_treasure/3)

    # Solve the system of equations for all variables
    solution = sympy.solve([eq1, eq2, eq3, eq4, eq5, eq6, eq7, eq8],
                           (p0, p1, p2, p3, p4, p5, p6, p7))

    # The result for the Main Entrance is in solution[p0]
    prob_main_entrance = solution[p0]
    
    # Extract the numerator and denominator for the final equation format
    numerator = prob_main_entrance.p
    denominator = prob_main_entrance.q

    print(f"The probability of reaching the Treasure Room from the Main Entrance is:")
    print(f"P(Main Entrance) = {numerator} / {denominator}")

solve_castle_escape()