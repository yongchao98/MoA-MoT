def solve_task():
    """
    This function identifies the true statements about probabilistic graphical models
    and prints them in the specified format.
    """
    # Statements are evaluated based on the principles of the Junction Tree algorithm.
    # A: False. Approximate inference methods can be more efficient.
    # B: False. The algorithm is not resilient to high treewidth; it's its main limitation.
    # C: True. The nodes are cliques, and the algorithm uses their joint distributions (potentials).
    # D: False. The premise is false.
    # E: True. High treewidth leads to intractably large clique tables.
    # F: False. Same as D.
    # G: True. The need to compute joint distributions over large cliques is why high treewidth is a problem.
    # H: False. The cost is exponential, not linear.
    # I: True. The cost grows exponentially with the size of the largest clique.
    # J: False. The cost changes dramatically.
    # L: True. This is a correct statement of the running intersection property.

    true_statements = ['C', 'E', 'G', 'I', 'L']

    # Format the output as a comma-separated list in curly brackets.
    output_string = "{" + ", ".join(true_statements) + "}"
    
    print(output_string)

solve_task()