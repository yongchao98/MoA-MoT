def solve_chow_groups():
    """
    This function computes the pairs (m(X), M(X)) for four given algebraic varieties.
    The values are derived from established theorems in algebraic geometry.
    """

    # For X_1, a genus 2 curve:
    # m(X_1) is the gonality of the curve, which is 2. This is edeg(w, C) for a Weierstrass point w.
    # M(X_1) is edeg(p, C) for a general point p, which is g+1 = 2+1 = 3.
    m_X1 = 2
    M_X1 = 3
    pair_X1 = (m_X1, M_X1)

    # For X_2, a general genus 7 curve:
    # m(X_2) is the gonality of a general genus 7 curve, which is floor((7+3)/2) = 5.
    # M(X_2) is edeg(p, C) for a general point p, which is g+1 = 7+1 = 8.
    m_X2 = 5
    M_X2 = 8
    pair_X2 = (m_X2, M_X2)

    # For X_3, an Enriques surface:
    # The Chow group of 0-cycles CH_0(X) is isomorphic to Z.
    # This implies [p] = [q] for any two points p, q.
    # Therefore, edeg(X, p) = 1 for all p.
    m_X3 = 1
    M_X3 = 1
    pair_X3 = (m_X3, M_X3)

    # For X_4, the Grassmannian G(3,6):
    # This is a rationally connected variety, so CH_0(X) is isomorphic to Z.
    # As with the Enriques surface, edeg(X, p) = 1 for all p.
    m_X4 = 1
    M_X4 = 1
    pair_X4 = (m_X4, M_X4)

    # Format the output string as requested
    result = f"({m_X1}, {M_X1}), ({m_X2}, {M_X2}), ({m_X3}, {M_X3}), ({m_X4}, {M_X4})"
    print(result)

solve_chow_groups()