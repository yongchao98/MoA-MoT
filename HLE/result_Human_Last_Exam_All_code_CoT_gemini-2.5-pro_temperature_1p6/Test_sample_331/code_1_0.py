def compute_edeg_pairs():
    """
    Computes the pairs (m(X), M(X)) for the four given varieties.

    The computation is based on established results in algebraic geometry.
    - edeg(X,p) is the smallest d > 0 such that for any q_1 in X,
      d[p] - [q_1] is rationally equivalent to an effective 0-cycle of degree d-1.
    - m(X) is the minimum of edeg(X,p) over all p in X.
    - M(X) is the supremum of edeg(X,p) over all p in X.
    """

    # For X_1, a genus 2 curve.
    # A genus 2 curve is hyperelliptic. For a hyperelliptic curve C of genus g,
    # edeg(C, p) = 2 for all points p.
    # So, m(X_1) = 2 and M(X_1) = 2.
    m_X1 = 2
    M_X1 = 2

    # For X_2, a general genus 7 curve.
    # For a general curve C of genus g, it is a known result that
    # edeg(C, p) = g for all points p.
    # Here, g=7. So, m(X_2) = 7 and M(X_2) = 7.
    m_X2 = 7
    M_X2 = 7

    # For X_3, an Enriques surface.
    # For an Enriques surface S, A_0(S) is Z/2Z. For any point p,
    # it is possible to find a point q such that [p]-[q] is the non-trivial element.
    # This implies that for d-1=1 (i.e., d=2), the condition for edeg can be satisfied.
    # This holds for all points, so m(X_3) = 2 and M(X_3) = 2.
    m_X3 = 2
    M_X3 = 2

    # For X_4, the Grassmannian G(3,6).
    # The Grassmannian is a rational variety. For any rational variety X,
    # the group of 0-cycles of degree 0, A_0(X), is trivial.
    # This means [p] is rationally equivalent to [q] for all p,q in X.
    # The condition for edeg(X,p) with d=1 is [p] = [q_1], which is true.
    # Thus, edeg(X,p) = 1 for all p. So m(X_4) = 1 and M(X_4) = 1.
    m_X4 = 1
    M_X4 = 1

    # Format the results as a string of pairs.
    result_string = f"({m_X1}, {M_X1}), ({m_X2}, {M_X2}), ({m_X3}, {M_X3}), ({m_X4}, {M_X4})"
    
    print(result_string)

compute_edeg_pairs()