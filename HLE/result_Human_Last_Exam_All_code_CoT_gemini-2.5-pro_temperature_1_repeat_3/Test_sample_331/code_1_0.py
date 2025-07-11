def compute_edeg_values():
    """
    This function computes the (m(X), M(X)) pairs for the four given varieties.
    The computation is based on established theorems in algebraic geometry concerning
    the Chow groups of these varieties.
    """

    # Case X_1: A genus 2 curve.
    # m(X) corresponds to Weierstrass points, M(X) to non-Weierstrass points.
    m_X1 = 2
    M_X1 = 3

    # Case X_2: A general genus 7 curve.
    # For a general curve of genus g, m(X) = M(X) = g+1. Here g=7.
    m_X2 = 8
    M_X2 = 8

    # Case X_3: An Enriques surface.
    # For an Enriques surface, CH_0(X)_0 = 0, so any two points are rationally equivalent.
    m_X3 = 1
    M_X3 = 1

    # Case X_4: A Grassmannian G(3,6).
    # As a rational variety, CH_0(X) is isomorphic to Z. Rational equivalence is equality of degree.
    m_X4 = 1
    M_X4 = 1

    # Format the result string as requested.
    # The requirement "output each number in the final equation" is interpreted as
    # constructing the string from the computed integer variables.
    result_string = f"({m_X1}, {M_X1}), ({m_X2}, {M_X2}), ({m_X3}, {M_X3}), ({m_X4}, {M_X4})"
    
    print(result_string)

compute_edeg_values()