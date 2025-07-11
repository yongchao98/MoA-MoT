def solve():
    """
    This function determines the complexity class for Luby's algorithm on three graph classes.
    """
    # Analysis for f_1(n) on a cycle of length n
    # The complexity is O(log n), and the lower bound is also Omega(log n)
    # because a cycle behaves like a long path.
    d1 = 9

    # Analysis for f_2(n) on any tree of degree at most 100
    # The complexity is O(log n). The worst-case tree is a path, which provides
    # an Omega(log n) lower bound.
    d2 = 9

    # Analysis for f_3(n) on any graph of degree at most 100
    # The complexity is known to be Theta(log n) for bounded-degree graphs.
    d3 = 9

    # The final answer is the concatenation of these digits.
    final_answer = f"{d1}{d2}{d3}"
    print(f"The analysis leads to the following complexity categories:")
    print(f"1) Cycles: f_1(n) = Theta(log n) -> Category {d1}")
    print(f"2) Trees (d_max<=100): f_2(n) = Theta(log n) -> Category {d2}")
    print(f"3) Graphs (d_max<=100): f_3(n) = Theta(log n) -> Category {d3}")
    print(f"\nThe resulting three-digit code is: {final_answer}")
    
    # The final output format required by the prompt
    print(f"<<<{final_answer}>>>")

solve()