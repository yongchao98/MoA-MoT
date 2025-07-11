def solve():
    """
    Identifies the true statements about the junction tree algorithm
    and prints them in the specified format.
    """
    # Based on the analysis, the following statements are true:
    # C: The junction tree uses the joint distributions within nodes (cliques).
    # E: High treewidth makes the junction tree algorithm impractical.
    # G: The reason for the impracticality with high treewidth is the use of
    #    joint distributions over large cliques.
    # I: The computational cost grows exponentially with the linear increase in the
    #    number of variables in the largest clique.
    # L: This is a correct description of the running intersection property.
    correct_statements = ['C', 'E', 'G', 'I', 'L']

    # Sort the list for consistent ordering.
    correct_statements.sort()

    # Format the output string as a comma-separated list within brackets.
    # The format specification is "{X, X, . . . }"
    output = "{" + ", ".join(correct_statements) + "}"

    print(output)

solve()