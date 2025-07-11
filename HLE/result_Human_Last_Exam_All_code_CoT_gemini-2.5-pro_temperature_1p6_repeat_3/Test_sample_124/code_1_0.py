def solve_task():
    """
    Analyzes the given statements about probabilistic graphical models and identifies the true ones.
    """

    # List of statements identified as true based on the analysis.
    # C: Correctly describes the use of joint distributions in clique nodes.
    # E: Correctly identifies high treewidth as a cause of impracticality.
    # G: Correctly links the lack of resilience to the use of joint distributions.
    # I: Correctly states the exponential relationship between cost and clique size.
    # L: Correctly describes the running intersection property.
    true_statements = ['C', 'E', 'G', 'I', 'L']

    # Sort the list for consistent ordering.
    true_statements.sort()

    # Format the output as a comma-separated list within brackets.
    # Example: {A, B, C}
    output_string = "{" + ", ".join(true_statements) + "}"

    print(output_string)

solve_task()