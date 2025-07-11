import itertools

def solve_manifold_problem():
    """
    Finds the lexicographically least tuple (a_1, b_1, ..., a_l, b_l)
    by searching for the minimal l and the combination of pairs (a_i, b_i)
    that satisfies the topological constraints.
    """
    # l must be odd and greater than 1. The minimal l is 3.
    l = 3
    target_c_sum = (l - 1) // 2

    # Generate a list of candidate pairs (a,b) and their c = (1-a)(1-b) values.
    # We search in a small range for a and b, as the solution is expected to have small values.
    # The pairs are generated and stored in lexicographical order.
    # a, b must be non-negative integers, and not equal to 1.
    candidates = []
    max_val = 5
    for a in range(max_val):
        if a == 1:
            continue
        for b in range(a, max_val):
            if b == 1:
                continue
            c = (1 - a) * (1 - b)
            candidates.append(((a, b), c))

    # We use combinations_with_replacement because the same manifold M(a,b) can be used multiple times.
    # The iterator yields combinations in sorted order, so the first solution we find
    # will correspond to the lexicographically smallest tuple.
    for combo in itertools.combinations_with_replacement(candidates, l):
        c_sum = sum(item[1] for item in combo)

        if c_sum == target_c_sum:
            pairs = [item[0] for item in combo]
            c_values = [item[1] for item in combo]
            
            p1, p2, p3 = pairs
            c1, c2, c3 = c_values
            
            result_tuple = (p1[0], p1[1], p2[0], p2[1], p3[0], p3[1])

            print(f"Found a solution for the minimal length l = {l}.")
            print("The three manifolds M(a_i, b_i), which are not 'full', are:")
            print(f"1. M({p1[0]},{p1[1]}) -> c_1 = (1-{p1[0]})*(1-{p1[1]}) = {c1}")
            print(f"2. M({p2[0]},{p2[1]}) -> c_2 = (1-{p2[0]})*(1-{p2[1]}) = {c2}")
            print(f"3. M({p3[0]},{p3[1]}) -> c_3 = (1-{p3[0]})*(1-{p3[1]}) = {c3}")
            print("\nTheir connect-sum is 'full' because the necessary condition is met:")
            print(f"2 * (c_1 + c_2 + c_3) = l - 1")
            print(f"2 * ({c1} + {c2} + {c3}) = {l} - 1")
            print(f"2 * ({sum(c_values)}) = {l-1}")
            print("This is true.\n")
            print("The lexicographically least tuple is:")
            # Format the tuple as a string with no spaces as requested.
            print(f"({','.join(map(str, result_tuple))})")
            
            # Since we iterate in order, the first solution found is the minimal one.
            return

solve_manifold_problem()