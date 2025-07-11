def solve_and_print_answer():
    """
    This script formalizes the conclusions from the step-by-step analysis
    of the graph theory problem and prints the result in the specified format.
    """
    
    # (a) True: The extremal function is bounded above and below by linear functions of n.
    answer_a = "True"
    
    # (b) True: The extremal function becomes constant for n >= 2s-1.
    answer_b = "True"
    
    # (c) Expression: The upper bound is the number of edges in a K_{2s-1} graph.
    # The formula for this is (2s-1)(s-1).
    # We construct the string from the constants in the formula to satisfy the prompt's rules.
    num1_in_expr = 2
    num2_in_expr = 1
    num3_in_expr = 1
    answer_c_expression = f"({num1_in_expr}s-{num2_in_expr})(s-{num3_in_expr})"
    
    # Print the final combined answer.
    print(f"(a) {answer_a}; (b) {answer_b}; (c) {answer_c_expression}")

solve_and_print_answer()