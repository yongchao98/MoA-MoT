def solve_chow_groups_problem():
    """
    This function encapsulates the theoretical results for the (m(X), M(X)) pairs.
    m(X) = min_{p in X} edeg(X,p)
    M(X) = sup_{p in X} edeg(X,p)
    """

    # X_1: genus 2 curve
    # A genus 2 curve is hyperelliptic.
    # edeg = 2 for Weierstrass points, edeg = 3 for other points.
    m_X1 = 2
    M_X1 = 3
    result_X1 = f"({m_X1}, {M_X1})"

    # X_2: general genus 7 curve
    # For a general curve of genus g, m(X) = M(X) = g + 1.
    # For g=7, this is 8.
    g = 7
    m_X2 = g + 1
    M_X2 = g + 1
    result_X2 = f"({m_X2}, {M_X2})"

    # X_3: Enriques surface
    # For an Enriques surface, the Chow group of 0-cycles of degree 0, A_0(X), is trivial.
    # This implies any two points are rationally equivalent, so edeg(X,p) = 1 for all p.
    m_X3 = 1
    M_X3 = 1
    result_X3 = f"({m_X3}, {M_X3})"

    # X_4: Grassmannian G(3,6)
    # G(3,6) is a rational variety. For rational varieties, A_0(X) is trivial.
    # This implies edeg(X,p) = 1 for all p.
    m_X4 = 1
    M_X4 = 1
    result_X4 = f"({m_X4}, {M_X4})"

    # Format the final output string
    final_answer = f"{result_X1}, {result_X2}, {result_X3}, {result_X4}"
    print(final_answer)

solve_chow_groups_problem()