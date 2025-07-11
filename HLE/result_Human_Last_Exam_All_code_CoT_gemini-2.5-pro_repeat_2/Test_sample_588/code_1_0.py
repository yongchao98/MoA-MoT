def solve_manifold_problem():
    """
    This function finds the lexicographically least tuple (a1, b1, ..., al, bl)
    satisfying the conditions of the problem.

    The problem simplifies to finding the lexicographically smallest tuple
    (a1, b1, a2, b2, a3, b3) such that:
    1. a_i, b_i are non-negative integers not equal to 1.
    2. (1-a1)(1-b1) + (1-a2)(1-b2) + (1-a3)(1-b3) = 1.
    """

    # The possible values for a_i and b_i are {0, 2, 3, 4, ...}
    # We only need to search a small range to find the minimal solution.
    search_range = [0] + list(range(2, 5))

    def get_k(a, b):
        """Calculates the term (1-a)(1-b) for the main equation."""
        return (1 - a) * (1 - b)

    # We search for the 6-tuple (a1, b1, a2, b2, a3, b3)
    # by iterating through possibilities in lexicographical order.
    # The first solution we find will be the smallest.
    for a1 in search_range:
        for b1 in search_range:
            k1 = get_k(a1, b1)
            for a2 in search_range:
                for b2 in search_range:
                    k2 = get_k(a2, b2)
                    for a3 in search_range:
                        for b3 in search_range:
                            k3 = get_k(a3, b3)
                            
                            # Check if the core equation is satisfied
                            if k1 + k2 + k3 == 1:
                                final_tuple = (a1, b1, a2, b2, a3, b3)
                                # Print the final tuple in the required format
                                print(f"({','.join(map(str, final_tuple))})")
                                return

solve_manifold_problem()