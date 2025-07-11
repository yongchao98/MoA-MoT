def find_true_statements():
    """
    This function identifies the correct statements about the Junction Tree algorithm
    and prints them in the specified format.
    
    A: False. Approximate inference can be more efficient.
    B: False. It is not resilient to high treewidth.
    C: True. Its nodes are cliques, and it uses their joint distributions (potentials).
    D: False. The premise is false.
    E: True. High treewidth makes it impractical.
    F: False. Duplicate of D.
    G: True. The need for joint distributions over large cliques is the reason for its intractability in high-treewidth graphs.
    H: False. The relationship is exponential, not linear.
    I: True. This correctly states the exponential complexity.
    J: False. Efficiency is highly dependent on the largest clique size.
    L: True. This is a correct description of the running intersection property.
    """
    
    # The list of letters corresponding to the true statements.
    true_statements = ['C', 'E', 'G', 'I', 'L']
    
    # Format the output as a comma-separated list enclosed in brackets.
    output = "{" + ", ".join(true_statements) + "}"
    
    print(output)

find_true_statements()
<<<{"C", "E", "G", "I", "L"}>>>