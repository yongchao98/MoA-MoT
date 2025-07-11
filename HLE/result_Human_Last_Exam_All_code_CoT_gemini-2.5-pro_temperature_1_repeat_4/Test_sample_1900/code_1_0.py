def solve_complexity_questions():
    """
    This function prints the computational complexity for the two questions based on graph theory principles.
    
    The problem described corresponds to the Hamiltonian Path Problem on a specific class of graphs.
    - Question A asks for the complexity of the decision problem.
    - Question B asks for the complexity of the search problem.
    
    The Hamiltonian Path Problem is NP-complete, even for the class of graphs described (locally connected).
    This means the best-known algorithms have exponential time complexity in terms of n, the number of red balls.
    We express this exponential complexity using big-O notation as O(e^n). The complexity is the same for both
    the decision and search versions of the problem.
    """
    
    # The complexity is exponential for both questions.
    # 'e' is Euler's number, a base for the natural logarithm. Using it in O(e^n) is a way
    # to denote exponential complexity.
    answer_a = "O(e^n)"
    answer_b = "O(e^n)"
    
    # The final equation is the combination of both answers separated by a semicolon.
    final_equation = f"${answer_a}$; ${answer_b}$"
    
    print(final_equation)

solve_complexity_questions()