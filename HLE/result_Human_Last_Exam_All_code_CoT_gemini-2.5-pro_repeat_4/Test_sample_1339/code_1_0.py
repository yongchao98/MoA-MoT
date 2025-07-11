def solve_group_theory_problem():
    """
    This function provides the solution to the group theory problem.
    """

    # Part (a): Existence and uniqueness of the minimal group.
    # This asks for the existence and uniqueness of the p-localization of G,
    # which is a standard result for this class of solvable groups.
    answer_a = "Yes"

    # Part (b): Maximum possible derived length of hat(G).
    # The derived length of G, dl(G), is at most n.
    # The derived length of the p-localization hat(G) is bounded by 2*dl(G) - 1.
    # The maximum possible value is therefore 2*n - 1.
    
    # The final equation is: max_derived_length = 2 * n - 1
    # We identify the numbers in this equation as requested.
    number_1_in_equation = 2
    number_2_in_equation = 1
    
    # We construct the expression string.
    expression_b = f"{number_1_in_equation}*n - {number_2_in_equation}"

    # Print the final answer in the specified format.
    print(f"(a) {answer_a}")
    print(f"(b) {expression_b}")

solve_group_theory_problem()