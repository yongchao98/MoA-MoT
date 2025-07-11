def solve_complexity_problem():
    """
    This function analyzes the computational complexity of the described game tasks
    and prints the answer in Big-O notation.
    """

    # The problem of finding a line that passes through all n red balls is equivalent
    # to the Hamiltonian Path Problem in a graph where balls are vertices and
    # adjacency represents neighboring balls of the same color.

    # Question A asks for the complexity of *deciding* if such a path exists.
    # The Hamiltonian Path Problem is a classic NP-complete problem. This means
    # there is no known algorithm that solves it in polynomial time. The best
    # known algorithms have a worst-case time complexity that is exponential
    # in the number of vertices, n. This is commonly expressed as O(2^n).
    complexity_A = "O(2^n)"

    # Question B asks for the complexity of *finding* the path, given it exists.
    # This is the search version of the problem. For NP-complete problems, the search
    # and decision versions are considered to be of similar difficulty. Algorithms
    # like the dynamic programming approach (based on Held-Karp) can solve both
    # the decision and search versions in O(n^2 * 2^n) time.
    # In terms of exponential complexity, this is also expressed as O(2^n).
    complexity_B = "O(2^n)"

    # The final answer is formatted as requested.
    print(f"{complexity_A}; {complexity_B}")

solve_complexity_problem()