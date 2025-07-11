def solve():
    """
    This function prints the computational complexity for the two questions.
    """
    # The complexity of deciding if a Hamiltonian path exists in the described graph.
    # The problem is NP-complete, and a common algorithm runs in exponential time.
    complexity_A = "O(n * 2^n)"

    # The complexity of finding the path, given one exists.
    # The same algorithm can be used to find the path with the same complexity.
    complexity_B = "O(n * 2^n)"

    # The final answer format is "A; B"
    print(f"{complexity_A}; {complexity_B}")

solve()