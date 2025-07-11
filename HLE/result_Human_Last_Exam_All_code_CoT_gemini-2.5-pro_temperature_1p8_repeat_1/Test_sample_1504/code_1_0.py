def solve_turran_extremal_problem():
    """
    Solves the Turán-type extremal problem parts (a), (b), and (c).
    The solution is derived from graph theory principles as explained above.
    The code constructs and prints the final answer string.
    """

    # Part (a) is True as ex(...) is O(n) from the max-degree constraint
    # and Omega(n) from a matching construction.
    answer_a = "True"

    # Part (b) is True as the number of edges is bounded by a constant
    # dependent only on s and t.
    answer_b = "True"

    # Part (c) requires an expression for the upper bound.
    # From e <= tau*Delta and tau <= 2*alpha', we get e <= 2(s-1)(t-1).
    # The numbers in the expression are derived as follows:
    s_coefficient = 2
    s_subtract_val = 1
    t_subtract_val = 1
    
    # We construct the expression string as requested.
    answer_c = f"{s_coefficient}*(s - {s_subtract_val})*(t - {t_subtract_val})"

    final_answer = f"(a) {answer_a}; (b) {answer_b}; (c) {answer_c}"
    print(final_answer)

solve_turran_extremal_problem()