def solve_puzzle():
    """
    Calculates the maximum possible number of digits in N based on the problem's rules.
    N uses at most 5 distinct digits.
    Let L(k) be the maximum length of a valid number with k distinct digits.
    """

    print("Step-by-step derivation of the maximum length L(k):")
    print("Let L(k) be the maximum possible number of digits in N using k distinct digits.")
    print("The recurrence relation is L(k) = 2 * L(k-1) + 1.")
    print("-" * 30)

    # Base case: k=1
    k = 1
    # For k=1, the only valid number is a single digit (e.g., '1'). '11' is invalid.
    l_k = 1
    print(f"For k = {k} digit:")
    print(f"L({k}) = {l_k}")
    print("-" * 30)

    # Iterate from k=2 to k=5
    for k in range(2, 6):
        l_prev = l_k
        l_k = 2 * l_prev + 1
        print(f"For k = {k} digits:")
        print(f"L({k}) = 2 * L({k-1}) + 1")
        print(f"L({k}) = 2 * {l_prev} + 1")
        print(f"L({k}) = {l_k}")
        print("-" * 30)

    print(f"The maximum possible number of digits in N using at most 5 distinct digits is L(5).")
    print(f"Final Answer: {l_k}")

solve_puzzle()