def solve_chow_group_problem():
    """
    This function stores and prints the computed values for (m(X), M(X))
    for the four given varieties.
    """

    # X_1: a genus 2 curve
    # For a genus 2 curve, m(X) = 3 and M(X) = 3.
    m_X1 = 3
    M_X1 = 3

    # X_2: a general genus 7 curve
    # For a general genus g curve, m(X) = M(X) = g+1. Here g=7.
    m_X2 = 8
    M_X2 = 8

    # X_3: an Enriques surface
    # For an Enriques surface, m(X) = 2 and M(X) = 2.
    m_X3 = 2
    M_X3 = 2

    # X_4: Grassmannian G(3,6)
    # For a rationally connected variety like G(3,6), m(X) = 1 and M(X) = 1.
    m_X4 = 1
    M_X4 = 1

    # The problem requires the output to be a sequence of pairs.
    # The f-string explicitly constructs the final string with each number.
    result_string = f"({m_X1}, {M_X1}), ({m_X2}, {M_X2}), ({m_X3}, {M_X3}), ({m_X4}, {M_X4})"

    print(result_string)

solve_chow_group_problem()