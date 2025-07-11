def solve_complexity_questions():
    """
    This function determines and prints the computational complexity for the two parts of the problem.
    """

    # Question A: What is the computational complexity of deciding whether I can draw a line through all the red balls in the game?
    #
    # The problem states that the graph of red balls is connected and locally connected.
    # A theorem in graph theory guarantees that such a graph (with n >= 2 vertices)
    # always contains a Hamiltonian path (a line through all balls).
    # For n=1, a path trivially exists.
    # Therefore, the decision is always "yes" and requires no computation on the graph itself.
    # This is a constant time operation.
    # The number in this complexity term is 1.
    complexity_A = "O(1)"

    # Question B: If I know that I can draw such a line, what is the computational complexity of finding a line that goes through all the red balls?
    #
    # This is the search version of the Hamiltonian Path problem. For the specific class of
    # connected, locally connected graphs, there are polynomial-time algorithms.
    # These algorithms are typically based on constructive proofs that iteratively extend a path.
    # The complexity of finding a Hamiltonian cycle (a closely related, slightly harder problem)
    # in such graphs has been shown to be O(n^3). This complexity bound applies to finding a path as well.
    # The number in this complexity term is 3.
    complexity_B = "O(n^3)"

    # The final answer should be in the format "O(ans_A); O(ans_B)"
    print(f"{complexity_A}; {complexity_B}")

solve_complexity_questions()