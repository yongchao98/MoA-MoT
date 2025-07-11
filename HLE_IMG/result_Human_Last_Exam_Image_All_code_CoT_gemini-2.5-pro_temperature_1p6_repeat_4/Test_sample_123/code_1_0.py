import sympy

def solve_bran_castle_escape():
    """
    Solves the Bran Castle Escape Challenge by setting up and solving a system of
    linear equations for the probabilities of reaching the Treasure Room.
    """
    # Define symbols for the probabilities of reaching the Treasure Room from each room.
    # The 'p' prefix stands for 'probability'.
    p_me, p_gh, p_ck, p_sp, p_l, p_kh, p_qr, p_tc = sympy.symbols(
        'p_me p_gh p_ck p_sp p_l p_kh p_qr p_tc'
    )

    # The probabilities for the absorbing states (the destinations).
    p_tr = 1  # Probability of reaching the treasure room from the treasure room is 1.
    p_vl = 0  # Probability of reaching the treasure room from the vampire's lair is 0.

    # Set up the system of linear equations based on the graph transitions.
    # For each room, the probability of winning from there is the sum of
    # (probability of taking an exit * probability of winning from the destination room).

    # 1. Main Entrance (ME) -> Great Hall (1/2), Castle Kitchen (1/2)
    eq1 = sympy.Eq(p_me, sympy.Rational(1, 2) * p_gh + sympy.Rational(1, 2) * p_ck)

    # 2. Great Hall (GH) -> Main Entrance (1/4), Secret Passage (1/4), Library (1/4), Knights' Hall (1/4)
    eq2 = sympy.Eq(p_gh, sympy.Rational(1, 4) * p_me + sympy.Rational(1, 4) * p_sp + sympy.Rational(1, 4) * p_l + sympy.Rational(1, 4) * p_kh)

    # 3. Castle Kitchen (CK) -> Main Entrance (1/2), Knights' Hall (1/2)
    eq3 = sympy.Eq(p_ck, sympy.Rational(1, 2) * p_me + sympy.Rational(1, 2) * p_kh)

    # 4. Secret Passage (SP) -> Great Hall (1/4), Torture Chamber (1/4), Treasure Room (1/4), Vampire's Lair (1/4)
    eq4 = sympy.Eq(p_sp, sympy.Rational(1, 4) * p_gh + sympy.Rational(1, 4) * p_tc + sympy.Rational(1, 4) * p_tr + sympy.Rational(1, 4) * p_vl)

    # 5. Library (L) -> Great Hall (1/3), Queen's Room (1/3), Knights' Hall (1/3)
    eq5 = sympy.Eq(p_l, sympy.Rational(1, 3) * p_gh + sympy.Rational(1, 3) * p_qr + sympy.Rational(1, 3) * p_kh)

    # 6. Knights' Hall (KH) -> Great Hall (1/4), Library (1/4), Queen's Room (1/4), Castle Kitchen (1/4)
    eq6 = sympy.Eq(p_kh, sympy.Rational(1, 4) * p_gh + sympy.Rational(1, 4) * p_l + sympy.Rational(1, 4) * p_qr + sympy.Rational(1, 4) * p_ck)

    # 7. Queen's Room (QR) -> Library (1/3), Knights' Hall (1/3), Treasure Room (1/3)
    eq7 = sympy.Eq(p_qr, sympy.Rational(1, 3) * p_l + sympy.Rational(1, 3) * p_kh + sympy.Rational(1, 3) * p_tr)

    # 8. Torture Chamber (TC) -> Secret Passage (1/2), Vampire's Lair (1/2)
    eq8 = sympy.Eq(p_tc, sympy.Rational(1, 2) * p_sp + sympy.Rational(1, 2) * p_vl)

    # Solve the system of 8 linear equations for the 8 unknown probabilities.
    solution = sympy.solve(
        [eq1, eq2, eq3, eq4, eq5, eq6, eq7, eq8],
        (p_me, p_gh, p_ck, p_sp, p_l, p_kh, p_qr, p_tc)
    )

    # Extract the probability for the Main Entrance, which is our starting point.
    result_fraction = solution[p_me]
    numerator = result_fraction.p
    denominator = result_fraction.q

    # Print the final result.
    print("The probability of eventually reaching the Treasure Room from the Main Entrance is given by the equation:")
    print(f"p_me = {numerator} / {denominator}")

solve_bran_castle_escape()