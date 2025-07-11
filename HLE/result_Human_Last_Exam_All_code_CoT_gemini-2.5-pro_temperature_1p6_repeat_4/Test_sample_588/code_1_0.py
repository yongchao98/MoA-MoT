import itertools

def solve_and_print():
    """
    This function implements the logic to find the lexicographically least tuple
    that solves the problem, based on the derived mathematical conditions.
    """
    # Based on the analysis, the problem reduces to finding the lexicographically
    # smallest tuple (a1, b1, a2, b2, a3, b3) that solves the equation:
    # (1-a1)(1-b1) + (1-a2)(1-b2) + (1-a3)(1-b3) = 1
    # under the constraint that a_i, b_i are non-negative integers not equal to 1.
    
    def pair_generator(limit):
        # Generates pairs (a, b) such that a,b >= 0, a,b != 1, a <= b,
        # ordered lexicographically.
        for a in range(limit):
            if a == 1:
                continue
            for b in range(a, limit):
                if b == 1:
                    continue
                yield (a, b)

    # We expect the smallest solution to involve small integers.
    search_limit = 5
    # Generate a list of possible pairs in lexicographical order.
    possible_pairs = list(pair_generator(search_limit))

    # We use combinations_with_replacement to select three pairs from our list.
    # This method ensures that the combinations are produced in an order that,
    # when concatenated, will yield lexicographically sorted tuples. Thus, the
    # first solution we find will be the minimal one.
    for p1, p2, p3 in itertools.combinations_with_replacement(possible_pairs, 3):
        x1 = (1 - p1[0]) * (1 - p1[1])
        x2 = (1 - p2[0]) * (1 - p2[1])
        x3 = (1 - p3[0]) * (1 - p3[1])

        if x1 + x2 + x3 == 1:
            # We have found the first, and therefore smallest, solution.
            final_tuple = p1 + p2 + p3

            # As requested, output the numbers in the final equation.
            print(f"The final equation is derived from l-1 = 2 * Sum[(1-a_i)(1-b_i)].")
            print(f"For minimal l=3, this is Sum[(1-a_i)(1-b_i)] = 1.")
            print(f"The equation with values is: ({x1}) + ({x2}) + ({x3}) = 1")
            
            # The final answer is the flat tuple, formatted as requested.
            print(str(final_tuple).replace(" ", ""))
            return

solve_and_print()