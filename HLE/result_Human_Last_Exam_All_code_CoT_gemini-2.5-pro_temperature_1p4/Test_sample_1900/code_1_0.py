def solve_complexity_questions():
    """
    This function prints the computational complexity for the two questions asked.
    
    The problem of deciding whether a line can be drawn through all n red balls
    is equivalent to the Hamiltonian Path Problem on the graph of red balls.
    The graph has n vertices and is connected and locally connected.

    Question A: Complexity of deciding if a Hamiltonian path exists.
    This problem is NP-complete. The complexity is exponential, O(e^n).

    Question B: Complexity of finding such a path.
    This is the search version of the same problem, which is NP-hard.
    The complexity is also exponential, O(e^n).
    """
    
    # The complexity for both Question A (decision) and Question B (search) is exponential.
    # We follow the format provided in the problem description.
    # "O(e^n)" represents exponential time complexity.
    answer_A = "O(e^n)"
    answer_B = "O(e^n)"
    
    # The final output is the two answers separated by a semicolon.
    final_answer = f"{answer_A}; {answer_B}"
    
    print(final_answer)

solve_complexity_questions()