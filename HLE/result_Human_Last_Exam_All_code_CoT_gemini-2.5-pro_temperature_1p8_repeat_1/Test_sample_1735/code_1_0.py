def solve():
    """
    Analyzes the statements and prints the letters of the true ones.
    """
    
    # Based on the analysis, the following statements are true:
    # A) The optimal solution can be found in O(n log n) time.
    # C) The optimal maximum cumulative impact is invariant under negation of all elements.
    # G) The optimal solution remains optimal if we append any number with magnitude larger than max{C(k)}.
    # H) There exists a sequence where the optimal solution requires alternating positive and negative numbers.
    
    true_statements = ['A', 'C', 'G', 'H']
    
    # The statements should be sorted in alphabetical order.
    # They are already sorted in the list above.
    
    print(", ".join(true_statements))

solve()
<<<A, C, G, H>>>