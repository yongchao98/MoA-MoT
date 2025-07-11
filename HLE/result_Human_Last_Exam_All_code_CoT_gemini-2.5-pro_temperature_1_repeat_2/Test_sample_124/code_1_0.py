def solve():
    """
    Analyzes the statements about the junction tree algorithm and identifies the correct ones.
    
    A: False. Approximate inference can be more efficient.
    B: False. The algorithm is not resilient to high treewidth.
    C: True. The nodes are cliques, and the algorithm uses their joint distributions.
    D: False. Based on a false premise.
    E: True. High treewidth makes the algorithm's complexity (which is exponential) impractical.
    F: False. Repetition of D.
    G: True. The need to handle joint distributions over large cliques (from high treewidth) is the reason for the lack of resilience.
    H: False. The complexity increases exponentially, not linearly.
    I: True. The computational cost is exponential in the size of the largest clique.
    J: False. The cost definitely changes.
    L: True. This is a correct statement of the running intersection property's implication.
    
    The correct statements are C, E, G, I, L.
    """
    correct_statements = ["C", "E", "G", "I", "L"]
    # The output format requires a comma-separated list enclosed in brackets.
    # The final answer format is specified as <<<answer content>>>.
    # The code will print the required format as a string.
    
    output = "{" + ", ".join(correct_statements) + "}"
    print(output)

solve()
print("<<<{" + ", ".join(["C", "E", "G", "I", "L"]) + "}>>>")