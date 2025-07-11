def solve_hss_multiplication():
    """
    Calculates the number of submatrix multiplications for an HSS matrix
    multiplication with a tree of a given depth.
    """
    depth = 4
    print(f"Calculating the number of submatrix multiplications for an HSS tree of depth {depth}.")
    print("The recurrence relation is M(d) = 2 * M(d-1) + 6, where d is the depth.\n")

    # Base case: For a tree of depth d=1 (a single leaf), there is 1 multiplication.
    m_prev = 1
    print(f"Base Case: For a tree of depth 1, the number of multiplications M(1) = {m_prev}")

    # Apply the recurrence relation for d = 2, 3, 4
    for d in range(2, depth + 1):
        m_curr = 2 * m_prev + 6
        print(f"Step d={d}: M({d}) = 2 * M({d-1}) + 6 = 2 * {m_prev} + 6 = {m_curr}")
        m_prev = m_curr

    final_answer = m_prev
    print(f"\nThus, for a Hierarchical Semi-separable tree with depth 4, a total of {final_answer} submatrix multiplications are performed.")

solve_hss_multiplication()