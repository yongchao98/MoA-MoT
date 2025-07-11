def solve_chow_group_problem():
    """
    This function computes and prints the pairs (m(X), M(X)) for four given varieties.
    The derivation of these values is based on principles of algebraic geometry.
    """

    # X_1: a genus 2 curve
    # m(X_1) = 2 (achieved at Weierstrass points)
    # M(X_1) = 3 (achieved at non-Weierstrass points)
    m1, M1 = 2, 3

    # X_2: a general genus 7 curve
    # m(X_2) = gonality - 1 = 4 - 1 = 3
    # M(X_2) = genus = 7
    m2, M2 = 3, 7

    # X_3: an Enriques surface
    # For Enriques surfaces, CH_0(X) is isomorphic to Z.
    # edeg(X, p) = 1 for all p.
    m3, M3 = 1, 1

    # X_4: G(3,6), a Grassmannian
    # Grassmannians are rationally connected, so CH_0(X) is isomorphic to Z.
    # edeg(X, p) = 1 for all p.
    m4, M4 = 1, 1

    # Format the result string as requested.
    # The format is (m(X_1), M(X_1)), (m(X_2), M(X_2)), (m(X_3), M(X_3)), (m(X_4), M(X_4))
    result_string = f"({m1}, {M1}), ({m2}, {M2}), ({m3}, {M3}), ({m4}, {M4})"

    print(result_string)

solve_chow_group_problem()