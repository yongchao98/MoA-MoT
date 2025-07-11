def solve_chow_groups():
    """
    This function computes the pairs (m(X), M(X)) for the four given varieties.
    The reasoning is based on standard results in algebraic geometry.

    1.  X_1: Genus 2 curve.
        edeg(p) is 2 if p is a Weierstrass point (6 of them), and 3 otherwise.
        m(X_1) = 2, M(X_1) = 3.

    2.  X_2: General genus 7 curve.
        For any point p on a general g=7 curve, edeg(p) = g+1 = 8.
        m(X_2) = 8, M(X_2) = 8.

    3.  X_3: Enriques surface.
        For any Enriques surface, 2([p]-[q])=0 in CH_0(S). This implies edeg(p) <= 2.
        Since CH_0(S)_0 is not always trivial, edeg(p) is not always 1.
        So m(X_3) = 2, M(X_3) = 2.

    4.  X_4: Grassmannian G(3,6).
        This is a rationally connected variety, so [p]=[q] for all p,q.
        This implies edeg(p) = 1 for all p.
        m(X_4) = 1, M(X_4) = 1.
    """

    # For X_1, a genus 2 curve
    m_X1 = 2
    M_X1 = 3

    # For X_2, a general genus 7 curve
    m_X2 = 8
    M_X2 = 8

    # For X_3, an Enriques surface
    m_X3 = 2
    M_X3 = 2

    # For X_4, the Grassmannian G(3,6)
    m_X4 = 1
    M_X4 = 1

    # Format the output string as requested
    result = f"({m_X1}, {M_X1}), ({m_X2}, {M_X2}), ({m_X3}, {M_X3}), ({m_X4}, {M_X4})"
    print(result)

solve_chow_groups()