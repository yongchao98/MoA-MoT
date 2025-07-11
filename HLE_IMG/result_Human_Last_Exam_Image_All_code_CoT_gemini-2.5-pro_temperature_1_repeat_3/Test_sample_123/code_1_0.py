import sympy

def solve_bran_castle_escape():
    """
    Solves the Bran Castle Escape Challenge by setting up and solving a system of linear equations.
    """
    # For clarity, let's map variable names to room names.
    # p0: Main Entrance
    # p1: Great Hall
    # p2: Castle Kitchen
    # p3: Secret Passage
    # p4: Library
    # p5: Knights' Hall
    # p6: Queen's Room
    # p7: Torture Chamber
    # The probability of success from the Treasure Room is 1.
    # The probability of success from the Vampire's Lair is 0.

    # 1. Define symbolic variables for the probabilities of reaching the Treasure Room from each room.
    p0, p1, p2, p3, p4, p5, p6, p7 = sympy.symbols('p0, p1, p2, p3, p4, p5, p6, p7')

    # 2. Set up the system of 8 linear equations based on the castle layout and exit probabilities.
    # Each equation represents: p_current_room = sum(prob_exit * p_next_room)

    # eq0: Main Entrance -> 1/2 Great Hall, 1/2 Castle Kitchen
    eq0 = sympy.Eq(p0, sympy.Rational(1, 2) * p1 + sympy.Rational(1, 2) * p2)

    # eq1: Great Hall -> 1/4 Main Entrance, 1/4 Secret Passage, 1/4 Library, 1/4 Knights' Hall
    eq1 = sympy.Eq(p1, sympy.Rational(1, 4) * p0 + sympy.Rational(1, 4) * p3 + sympy.Rational(1, 4) * p4 + sympy.Rational(1, 4) * p5)

    # eq2: Castle Kitchen -> 1/2 Main Entrance, 1/2 Knights' Hall
    eq2 = sympy.Eq(p2, sympy.Rational(1, 2) * p0 + sympy.Rational(1, 2) * p5)

    # eq3: Secret Passage -> 1/4 Great Hall, 1/4 Torture Chamber, 1/4 Treasure Room (prob=1), 1/4 Vampire's Lair (prob=0)
    eq3 = sympy.Eq(p3, sympy.Rational(1, 4) * p1 + sympy.Rational(1, 4) * p7 + sympy.Rational(1, 4) * 1)

    # eq4: Library -> 1/3 Great Hall, 1/3 Queen's Room, 1/3 Knights' Hall
    eq4 = sympy.Eq(p4, sympy.Rational(1, 3) * p1 + sympy.Rational(1, 3) * p6 + sympy.Rational(1, 3) * p5)

    # eq5: Knights' Hall -> 1/4 Great Hall, 1/4 Library, 1/4 Queen's Room, 1/4 Castle Kitchen
    eq5 = sympy.Eq(p5, sympy.Rational(1, 4) * p1 + sympy.Rational(1, 4) * p4 + sympy.Rational(1, 4) * p6 + sympy.Rational(1, 4) * p2)

    # eq6: Queen's Room -> 1/3 Library, 1/3 Knights' Hall, 1/3 Treasure Room (prob=1)
    eq6 = sympy.Eq(p6, sympy.Rational(1, 3) * p4 + sympy.Rational(1, 3) * p5 + sympy.Rational(1, 3) * 1)

    # eq7: Torture Chamber -> 1/2 Secret Passage, 1/2 Vampire's Lair (prob=0)
    eq7 = sympy.Eq(p7, sympy.Rational(1, 2) * p3)

    # 3. Solve the system of equations.
    equations = [eq0, eq1, eq2, eq3, eq4, eq5, eq6, eq7]
    solution = sympy.solve(equations, (p0, p1, p2, p3, p4, p5, p6, p7))

    # 4. Extract and print the probability for the Main Entrance (p0).
    prob_main_entrance = solution[p0]
    numerator = prob_main_entrance.p
    denominator = prob_main_entrance.q

    print("The system of linear equations being solved is:")
    for i, eq in enumerate(equations):
        print(f"  Eq{i}: {eq}")
    
    print("\nThe probability of eventually reaching the Treasure Room from the Main Entrance is:")
    print(f"P(Main Entrance) = {numerator}/{denominator}")
    
    # As requested: "output each number in the final equation!"
    # The final equation is P = numerator / denominator. The numbers are the numerator and the denominator.
    print(f"\nFinal Equation: P = {numerator} / {denominator}")

solve_bran_castle_escape()