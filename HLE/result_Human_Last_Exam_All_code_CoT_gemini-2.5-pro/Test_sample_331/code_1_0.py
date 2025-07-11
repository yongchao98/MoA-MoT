def compute_edeg_pairs():
    """
    This function computes and prints the pairs (m(X), M(X)) for the four given varieties.
    The values are pre-calculated based on principles of algebraic geometry.
    """

    # X_1: a genus 2 curve
    # m(X_1) = 2 (for Weierstrass points)
    # M(X_1) = 3 (for non-Weierstrass points)
    m_X1 = 2
    M_X1 = 3

    # X_2: a general genus 7 curve
    # For a general curve of genus g, edeg(p) = g+1 for all p.
    # g = 7, so edeg(p) = 8.
    m_X2 = 8
    M_X2 = 8

    # X_3: an Enriques surface
    # For an Enriques surface S, CH_0(S)_0 = 0, so [p]=[q] for all p,q.
    # This implies edeg(p) = 1 for all p.
    m_X3 = 1
    M_X3 = 1

    # X_4: G(3,6), a Grassmannian
    # As a rational variety, CH_0(G(3,6))_0 = 0.
    # The reasoning is the same as for the Enriques surface.
    m_X4 = 1
    M_X4 = 1

    # Print the formatted output string
    print(f"({m_X1}, {M_X1}), ({m_X2}, {M_X2}), ({m_X3}, {M_X3}), ({m_X4}, {M_X4})")

compute_edeg_pairs()