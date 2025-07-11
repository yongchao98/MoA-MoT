def solve_graph_theory_problem():
    """
    This function provides the solution to the Turán-type extremal problem.
    The analysis for each part is as follows:

    (a) True: If G is not a union of K2's, it has a component with >=3 vertices.
    In a G-free, K_{1,t}-ind-free graph, the max degree is bounded by a constant
    related to a Ramsey number, making the number of edges linear in n, i.e., Θ(n).

    (b) True: An sK2-free graph has a matching number at most s-1. This implies
    all edges are incident to a set of at most 2(s-1) vertices. The K_{1,t}-ind
    condition further limits the number of edges, making the total edge count
    bounded by a constant independent of n.

    (c) The tight upper bound is generally given by the construction K_{s-1} ∨ I_{t-1}.
    The number of edges in this graph is (s-1 choose 2) + (s-1)(t-1).
    This expression is formatted as a string below.
    """

    answer_a = "True"
    answer_b = "True"

    # The formula for part (c) is (s-1)*(s-2)/2 + (s-1)*(t-1).
    # The prompt requires that each number in the equation be outputted.
    # To satisfy this, the string is constructed from components showing the numbers.
    num_1_str = "1"
    num_2_str = "2"

    term1_numerator = f"(s - {num_1_str}) * (s - {num_2_str})"
    term1 = f"{term1_numerator} / {num_2_str}"

    term2_factor1 = f"(s - {num_1_str})"
    term2_factor2 = f"(t - {num_1_str})"
    term2 = f"{term2_factor1} * {term2_factor2}"

    answer_c = f"{term1} + {term2}"

    final_answer_string = f"(a) {answer_a}; (b) {answer_b}; (c) {answer_c}"
    print(final_answer_string)


solve_graph_theory_problem()