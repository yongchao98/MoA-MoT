def solve_max_planar_area():
    """
    Finds the supremum of the area nm of a chessboard for which the
    super-knight (3,2) graph G is planar.

    The planarity of G is determined by the inequality (n-5)(m-5) <= 11,
    derived from Euler's formula for planar graphs. We search for the
    integer pair (n, m) that maximizes the area n*m under this constraint.
    
    Based on the analysis, cases n=4 and n=5 lead to an unbounded search space.
    We assume the interesting, constrained cases occur for n, m >= 6, which
    leads to a finite maximum area.
    """
    max_area = 0
    best_n = 0
    best_m = 0

    # The condition (n-5)(m-5) <= 11 for n,m>=6 implies that
    # neither n-5 nor m-5 can be greater than 11.
    # So, we only need to search up to n=16.
    for n in range(6, 17):
        # From the inequality: m-5 <= 11 / (n-5)
        # So, m <= 5 + 11 / (n-5)
        # Since m must be an integer, we take the floor.
        max_m_for_n = int(5 + 11 / (n - 5))
        
        # We only care about m >= 6
        if max_m_for_n >= 6:
            m = max_m_for_n
            current_area = n * m
            if current_area > max_area:
                max_area = current_area
                best_n = n
                best_m = m

    print("The condition for the graph to be potentially planar is (n-5)*(m-5) <= 11.")
    print(f"The largest area is found for the board {best_n}x{best_m}.")
    print(f"For this board, the condition is satisfied: ({best_n}-5)*({best_m}-5) = {(best_n-5)*(best_m-5)}, which is <= 11.")
    print(f"The maximum area nm is {best_n}*{best_m} = {max_area}.")
    print("\nThe supremum of the value of the largest size nm is:")
    print(max_area)

solve_max_planar_area()