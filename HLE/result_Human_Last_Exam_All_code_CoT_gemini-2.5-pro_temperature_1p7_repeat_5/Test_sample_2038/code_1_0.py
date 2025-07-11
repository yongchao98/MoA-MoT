def solve_knot_problem():
    """
    This function finds the number of 2-bridge knots with crossing number at most 13
    that admit two disjoint non-parallel minimal genus Seifert surfaces.

    Based on knot theory, such knots are precisely the T(2, q) torus knots, where q is an odd integer
    with |q| > 1. The crossing number of T(2, q) is |q|.
    Since knots and their mirror images are considered the same, we only need to count positive odd integers q
    such that 1 < q <= 13.
    """
    max_crossing_number = 13
    # The condition is that q is an odd integer with |q| > 1.
    # We consider q > 0 due to mirror images being non-distinct.
    # So we start our search from q = 3.
    lower_bound = 3

    knot_q_values = []
    for q in range(lower_bound, max_crossing_number + 1):
        # We are looking for odd integers q.
        if q % 2 != 0:
            knot_q_values.append(q)

    count = len(knot_q_values)

    print(f"The knots satisfying the condition are the T(2, q) torus knots where q is an odd integer.")
    print(f"The crossing number is q, which must be at most {max_crossing_number}.")
    print(f"The possible values for q are: {', '.join(map(str, knot_q_values))}.")
    print("\nTo find the total count, we sum 1 for each valid knot found.")

    # Building and printing the final equation as requested.
    equation_parts = ["1"] * count
    equation_str = " + ".join(equation_parts)
    print(f"Final Equation: {equation_str} = {count}")


solve_knot_problem()