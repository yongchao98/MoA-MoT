import sys

def solve():
    """
    This function prints the computational complexity for the two questions.
    
    The problem described is equivalent to finding a Hamiltonian path in a specific type of graph.
    - Question A asks for the complexity of the decision problem (does a path exist?).
    - Question B asks for the complexity of the search problem (find the path).

    Both the decision and search versions of the Hamiltonian Path problem are NP-complete/NP-hard,
    even for the class of graphs described. The best-known algorithms run in exponential time.
    A standard algorithm for this problem has a complexity of O(n^2 * 2^n).
    """

    # The complexity is O(n^power_of_n * base^n)
    power_of_n = 2
    base_of_exponent = 2

    # The complexity is the same for both the decision problem (A) and the search problem (B).
    complexity_string = f"$O(n^{power_of_n} \\cdot {base_of_exponent}^n)$"

    # Print the answer in the requested format "A; B"
    print(f"{complexity_string}; {complexity_string}")

solve()