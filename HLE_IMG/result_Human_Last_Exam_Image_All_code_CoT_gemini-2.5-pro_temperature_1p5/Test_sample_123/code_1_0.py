import sympy

def solve_bran_castle_escape():
    """
    Solves the Bran Castle escape probability problem using a system of linear equations.
    """
    # Define the unknown probabilities as symbolic variables
    p_me, p_gh, p_ck, p_sp, p_l, p_kh, p_qr, p_tc = sympy.symbols(
        'p_me p_gh p_ck p_sp p_l p_kh p_qr p_tc'
    )

    # Known probabilities for absorbing states
    p_tr = 1  # Treasure Room
    p_vl = 0   # Vampire's Lair

    # Define the system of linear equations based on the castle layout
    # We rearrange each equation to the form f(vars) = 0
    # The numbers in the equations are kept as fractions to ensure precision
    eq1 = sympy.Eq(p_me, sympy.S(1)/2 * p_gh + sympy.S(1)/2 * p_ck)
    eq2 = sympy.Eq(p_gh, sympy.S(1)/4 * p_me + sympy.S(1)/4 * p_sp + sympy.S(1)/4 * p_l + sympy.S(1)/4 * p_kh)
    eq3 = sympy.Eq(p_ck, sympy.S(1)/2 * p_me + sympy.S(1)/2 * p_kh)
    eq4 = sympy.Eq(p_sp, sympy.S(1)/4 * p_gh + sympy.S(1)/4 * p_tc + sympy.S(1)/4 * p_tr + sympy.S(1)/4 * p_vl)
    eq5 = sympy.Eq(p_l, sympy.S(1)/3 * p_gh + sympy.S(1)/3 * p_qr + sympy.S(1)/3 * p_kh)
    eq6 = sympy.Eq(p_kh, sympy.S(1)/4 * p_gh + sympy.S(1)/4 * p_l + sympy.S(1)/4 * p_qr + sympy.S(1)/4 * p_ck)
    eq7 = sympy.Eq(p_qr, sympy.S(1)/3 * p_l + sympy.S(1)/3 * p_kh + sympy.S(1)/3 * p_tr)
    eq8 = sympy.Eq(p_tc, sympy.S(1)/2 * p_sp + sympy.S(1)/2 * p_vl)

    # Print the system of equations being solved
    print("System of Equations to be Solved:")
    equations_text = [
        "p_ME = 1/2 * p_GH + 1/2 * p_CK",
        "p_GH = 1/4 * p_ME + 1/4 * p_SP + 1/4 * p_L + 1/4 * p_KH",
        "p_CK = 1/2 * p_ME + 1/2 * p_KH",
        "p_SP = 1/4 * p_GH + 1/4 * p_TC + 1/4",
        "p_L  = 1/3 * p_GH + 1/3 * p_QR + 1/3 * p_KH",
        "p_KH = 1/4 * p_GH + 1/4 * p_L + 1/4 * p_QR + 1/4 * p_CK",
        "p_QR = 1/3 * p_L + 1/3 * p_KH + 1/3",
        "p_TC = 1/2 * p_SP"
    ]
    for eq_text in equations_text:
        print(eq_text)
    print("-" * 30)

    # Solve the system of equations
    solution = sympy.solve(
        [eq1, eq2, eq3, eq4, eq5, eq6, eq7, eq8],
        [p_me, p_gh, p_ck, p_sp, p_l, p_kh, p_qr, p_tc]
    )

    # Extract the probability for the Main Entrance
    prob_main_entrance = solution[p_me]
    numerator, denominator = prob_main_entrance.p, prob_main_entrance.q

    # Print the final result
    print(f"The probability of reaching the Treasure Room from the Main Entrance is an irreducible fraction.")
    print(f"Final Equation: P(Main Entrance) = {numerator} / {denominator}")
    print(f"The numerator of the final equation is: {numerator}")
    print(f"The denominator of the final equation is: {denominator}")
    print(f"The final probability is: {prob_main_entrance}")

if __name__ == '__main__':
    solve_bran_castle_escape()