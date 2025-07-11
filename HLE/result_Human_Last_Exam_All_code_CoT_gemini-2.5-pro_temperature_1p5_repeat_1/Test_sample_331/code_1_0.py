def solve_chow_groups():
    """
    Computes and prints the (m(X), M(X)) pairs for the four given varieties.
    """
    # X_1: genus 2 curve
    m_X1 = 2
    M_X1 = 3

    # X_2: general genus 7 curve
    m_X2 = 5
    M_X2 = 5

    # X_3: Enriques surface
    m_X3 = 2
    M_X3 = "infinity"

    # X_4: Grassmannian G(3,6)
    m_X4 = 1
    M_X4 = 1
    
    # Format the final output string as requested.
    result = f"({m_X1}, {M_X1}), ({m_X2}, {M_X2}), ({m_X3}, {M_X3}), ({m_X4}, {M_X4})"
    
    print(result)

solve_chow_groups()