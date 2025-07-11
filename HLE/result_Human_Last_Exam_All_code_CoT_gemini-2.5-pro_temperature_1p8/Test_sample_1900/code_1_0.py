def solve_complexity():
    """
    This function provides the computational complexity for the two questions
    based on the analysis of the graph problem.
    """
    # The complexity is determined by the Hamiltonian Path problem for n vertices.
    # The best known algorithms have a complexity of O(n^2 * 2^n).
    complexity = "O(n^2 * 2^n)"
    
    # Question A asks for the complexity of deciding if a path exists.
    complexity_A = complexity
    
    # Question B asks for the complexity of finding the path.
    complexity_B = complexity
    
    # The problem asks to separate the answers with a semicolon.
    print(f"{complexity_A}; {complexity_B}")

solve_complexity()