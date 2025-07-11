def solve_complexity_questions():
    """
    This function prints the computational complexity for the two questions based on the analysis.
    """
    # The complexity for Question A: Deciding if a path exists.
    # As explained, the properties of the graph guarantee a path exists.
    # Thus, the decision is trivial and can be done in constant time.
    complexity_A = "O(1)"

    # The complexity for Question B: Finding the path.
    # In the absence of a known polynomial-time algorithm for this specific type of graph,
    # we rely on general algorithms for the Hamiltonian Path problem. The Held-Karp
    # algorithm is a standard choice with exponential complexity. The numbers in the
    # expression are 2 (for the base of the exponent) and 2 (for the exponent of n).
    complexity_B = "O(n^2 * 2^n)"

    # The final answer format is the two complexities separated by a semicolon.
    final_answer = f"{complexity_A}; {complexity_B}"
    
    print(final_answer)

solve_complexity_questions()