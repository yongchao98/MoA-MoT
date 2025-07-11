def solve_feynman_diagram_count():
    """
    This function finds a(3), the number of non-vanishing Feynman diagrams of order 2*3=6.
    The values for a(n) are from a known sequence (OEIS A000599).
    """
    # a(n) is the number of diagrams of order 2n.
    # The sequence is known for the first few terms:
    # a(1) = 1
    # a(2) = 2
    # a(3) = 10
    # a(4) = 74
    # a(5) = 706
    feynman_sequence = {
        1: 1,
        2: 2,
        3: 10,
        4: 74,
        5: 706
    }

    n = 3
    
    if n in feynman_sequence:
        result = feynman_sequence[n]
        # The final equation is a(3) = 10
        # Outputting each number in the final equation as requested.
        print(f"a({n}) = {result}")
    else:
        print(f"The value for a({n}) is not stored.")

solve_feynman_diagram_count()