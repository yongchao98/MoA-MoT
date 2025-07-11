import sympy

def solve_bran_castle_escape():
    """
    Solves the Bran Castle Escape Challenge by setting up and solving a system of linear equations.
    """
    # Define the 8 unknown probabilities as symbols
    # Each symbol represents the probability of reaching the Treasure Room from that location.
    p_me, p_gh, p_sp, p_l, p_tc, p_qr, p_kh, p_ck = sympy.symbols(
        'p_me p_gh p_sp p_l p_tc p_qr p_kh p_ck'
    )

    # Probabilities for the absorbing states (the end points)
    p_tr = 1  # Probability of success from the Treasure Room is 1
    p_vl = 0  # Probability of success from the Vampire's Lair is 0

    # We create a system of 8 linear equations based on the transition probabilities
    # from each room. The structure of each equation is:
    # p_current_room = sum(probability_to_next_room * p_next_room)

    # Equation 1: Main Entrance (Exits to Great Hall, Castle Kitchen)
    eq1 = sympy.Eq(p_me, sympy.Rational(1, 2) * p_gh + sympy.Rational(1, 2) * p_ck)

    # Equation 2: Great Hall (Exits to Main Entrance, Secret Passage, Library, Knights' Hall)
    eq2 = sympy.Eq(p_gh, sympy.Rational(1, 4) * p_me + sympy.Rational(1, 4) * p_sp + sympy.Rational(1, 4) * p_l + sympy.Rational(1, 4) * p_kh)

    # Equation 3: Secret Passage (Exits to Great Hall, Torture Chamber, Treasure Room, Vampire's Lair)
    eq3 = sympy.Eq(p_sp, sympy.Rational(1, 4) * p_gh + sympy.Rational(1, 4) * p_tc + sympy.Rational(1, 4) * p_tr + sympy.Rational(1, 4) * p_vl)

    # Equation 4: Library (Exits to Great Hall, Queen's Room, Knights' Hall)
    eq4 = sympy.Eq(p_l, sympy.Rational(1, 3) * p_gh + sympy.Rational(1, 3) * p_qr + sympy.Rational(1, 3) * p_kh)

    # Equation 5: Torture Chamber (Exits to Secret Passage, Vampire's Lair)
    eq5 = sympy.Eq(p_tc, sympy.Rational(1, 2) * p_sp + sympy.Rational(1, 2) * p_vl)

    # Equation 6: Queen's Room (Exits to Library, Knights' Hall, Treasure Room)
    eq6 = sympy.Eq(p_qr, sympy.Rational(1, 3) * p_l + sympy.Rational(1, 3) * p_kh + sympy.Rational(1, 3) * p_tr)

    # Equation 7: Knights' Hall (Exits to Great Hall, Library, Queen's Room, Castle Kitchen)
    eq7 = sympy.Eq(p_kh, sympy.Rational(1, 4) * p_gh + sympy.Rational(1, 4) * p_l + sympy.Rational(1, 4) * p_qr + sympy.Rational(1, 4) * p_ck)

    # Equation 8: Castle Kitchen (Exits to Main Entrance, Knights' Hall)
    eq8 = sympy.Eq(p_ck, sympy.Rational(1, 2) * p_me + sympy.Rational(1, 2) * p_kh)

    # Define the list of variables we want to solve for
    variables = [p_me, p_gh, p_sp, p_l, p_tc, p_qr, p_kh, p_ck]

    # Define the list of equations to solve
    equations = [eq1, eq2, eq3, eq4, eq5, eq6, eq7, eq8]

    # Use sympy's linear equation solver to find the values of the variables
    solution = sympy.linsolve(equations, variables)

    # The solver returns a set containing one tuple of solutions.
    # We extract this tuple.
    solution_tuple = list(solution)[0]

    # Find the probability for the Main Entrance from the solution tuple
    prob_main_entrance = solution_tuple[variables.index(p_me)]

    # Extract the numerator and denominator for the final answer
    numerator = prob_main_entrance.p
    denominator = prob_main_entrance.q

    print(f"The exact probability of reaching the Treasure Room from the Main Entrance is an irreducible fraction.")
    print(f"Final Equation: P(Success) = {numerator} / {denominator}")
    print(f"The numerator is: {numerator}")
    print(f"The denominator is: {denominator}")

if __name__ == '__main__':
    solve_bran_castle_escape()