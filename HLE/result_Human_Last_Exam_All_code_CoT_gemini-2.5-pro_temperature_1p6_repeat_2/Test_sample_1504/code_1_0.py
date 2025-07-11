def solve_extremal_problem():
    """
    This function formulates and prints the solution to the Turan-type extremal problem.
    """
    # Part (a): Based on Ramsey's theorem and the properties of the forbidden graphs,
    # the extremal function grows linearly with n.
    answer_a = "True"

    # Part (b): Based on properties of matchings (Gallai's lemma) and the induced subgraph
    # condition, the number of edges is bounded by a constant independent of n.
    answer_b = "True"

    # Part (c): The upper bound is derived from the proof for part (b).
    # The bound is binom(2(s-1), 2) + 2(s-1)(t-1), which simplifies to (s-1)(2s+2t-5).
    # We will represent this expression as a formatted string.
    s = 's'
    t = 't'
    expression_c = f"({s} - 1) * (2*{s} + 2*{t} - 5)"

    # Print the final answer in the required format.
    # The final equation has its terms determined by the symbolic variables s and t.
    print(f"(a) {answer_a}; (b) {answer_b}; (c) {expression_c}")

# Execute the function to get the answer.
solve_extremal_problem()