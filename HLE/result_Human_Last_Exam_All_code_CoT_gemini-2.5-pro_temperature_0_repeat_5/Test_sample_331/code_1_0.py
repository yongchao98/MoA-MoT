def compute_edeg_pairs():
    """
    Computes the pairs (m(X), M(X)) for the four given varieties.

    The computation is based on established theorems in algebraic geometry.
    - X_1 (genus 2 curve): For any genus 2 curve, edeg(X,p) is always 2.
    - X_2 (general genus 7 curve): For a general curve of genus g, edeg(X,p) = g for all points p.
    - X_3 (Enriques surface): For an Enriques surface, the map from the surface to its
      group of 0-cycles of degree 0 is surjective, which implies edeg(X,p) = 2 for all p.
    - X_4 (Grassmannian G(3,6)): As a rational variety, all points are rationally
      equivalent, which implies edeg(X,p) = 1 for all p.
    """

    # For X_1, a genus 2 curve
    m_X1 = 2
    M_X1 = 2

    # For X_2, a general genus 7 curve
    m_X2 = 7
    M_X2 = 7

    # For X_3, an Enriques surface
    m_X3 = 2
    M_X3 = 2

    # For X_4, the Grassmannian G(3,6)
    m_X4 = 1
    M_X4 = 1

    # Format the output string as requested
    result = f"({m_X1}, {M_X1}), ({m_X2}, {M_X2}), ({m_X3}, {M_X3}), ({m_X4}, {M_X4})"
    print(result)

if __name__ == "__main__":
    compute_edeg_pairs()