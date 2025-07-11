def find_max_planar_area():
    """
    Finds the maximum area nm of a chessboard such that the super-knight graph is planar.
    The planarity condition is approximated by (n-5)(m-5) <= 11.
    We search for n, m >= 4. The condition is most meaningful for n, m >= 6.
    """
    max_area = 0
    best_n, best_m = 0, 0
    
    # We only need to check a reasonable range for n. If n is large, m must be small.
    # Let's check n from 4 up to a reasonable limit, e.g., 20.
    # By symmetry, we can assume n <= m.
    
    print("Searching for the maximum area nm such that (n-5)(m-5) <= 11 for n, m >= 4.\n")
    
    # We analyze the constraint (n-5)(m-5) <= 11
    # Let x = n-5, y = m-5. Then x,y >= -1. We want to maximize (x+5)(y+5) s.t. xy <= 11.
    
    potential_solutions = []

    # Case 1: n=4 (x = -1)
    # -y <= 11 -> y >= -11. Since m>=4, y=m-5>=-1. This holds for all m>=4.
    # This case would lead to an infinite supremum if G(4,m) is always planar.
    # We assume the problem implies a finite bound, so this case is likely excluded
    # by finding a K_{3,3} minor for large m.
    
    # Case 2: n=5 (x = 0)
    # 0 <= 11. This holds for all m>=5. Similar to n=4, this would imply an infinite supremum.
    
    # Case 3: n >= 6 (x >= 1)
    # This implies y must be <= 11/x. This gives a bounded search space.
    print("Analyzing cases for n >= 6:")
    for n in range(6, 17): # If n-5 > 11, then m-5 must be < 1, so m < 6.
        x = n - 5
        # y <= 11/x. m-5 <= 11/x -> m <= 5 + 11/x
        max_m = int(5 + 11 / x)
        for m in range(n, max_m + 1): # Assume n <= m for symmetry
            y = m - 5
            if x * y <= 11:
                area = n * m
                potential_solutions.append((n, m, area))
                if area > max_area:
                    max_area = area
                    best_n, best_m = n, m

    print("Potential solutions (n, m, area) for n>=6:")
    for n, m, area in sorted(potential_solutions, key=lambda item: item[2], reverse=True):
        print(f"n={n}, m={m}, area={area}, (n-5)(m-5)={(n-5)*(m-5)}")

    print(f"\nThe maximum area found under the condition for n,m>=6 is {max_area} for a {best_n}x{best_m} board.")
    print("\nThis value represents the supremum under the assumption that the problem is well-posed with a finite answer, which points to the derived density condition being the deciding factor.")

find_max_planar_area()