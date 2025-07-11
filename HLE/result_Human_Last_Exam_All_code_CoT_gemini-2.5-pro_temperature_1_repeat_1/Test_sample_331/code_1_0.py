def solve_variety_parameters():
    """
    This function computes and prints the (m(X), M(X)) pairs for the four given varieties.
    The values are derived from known results in algebraic geometry.
    """

    # For X_1, a genus 2 curve, m(X_1) = 2 and M(X_1) = 2.
    # A genus 2 curve is hyperelliptic, and it can be shown that for any point p,
    # edeg(X_1, p) = 2.
    m_X1 = 2
    M_X1 = 2
    pair_X1 = (m_X1, M_X1)

    # For X_2, a general genus 7 curve, m(X_2) = 3 and M(X_2) = 7.
    # M(X_2) is the genus g=7.
    # m(X_2) = g - gon(C) + 1 = 7 - 5 + 1 = 3.
    m_X2 = 3
    M_X2 = 7
    pair_X2 = (m_X2, M_X2)

    # For X_3, an Enriques surface, m(X_3) = 2 and M(X_3) = 2.
    # The Chow group of 0-cycles of degree 0, A_0(X_3), is a non-trivial finite group,
    # and the map from the surface to this group is surjective. This implies edeg(X_3, p) = 2 for all p.
    m_X3 = 2
    M_X3 = 2
    pair_X3 = (m_X3, M_X3)

    # For X_4, the Grassmannian G(3,6), m(X_4) = 1 and M(X_4) = 1.
    # G(3,6) is a rational variety, for which A_0(X_4) is trivial.
    # This means [p] = [q] for all points p,q, which makes edeg(X_4, p) = 1.
    m_X4 = 1
    M_X4 = 1
    pair_X4 = (m_X4, M_X4)

    # The problem asks to format the answer as a sequence of pairs.
    # The f-string below constructs the output with each number from the calculations above.
    result_string = f"{pair_X1}, {pair_X2}, {pair_X3}, {pair_X4}"
    
    print(result_string)

solve_variety_parameters()