def solve():
    """
    This function determines the category for each of the three scenarios based on the analysis.
    The analysis concludes that for each case, a path graph P_n serves as a worst-case example,
    leading to a number of steps proportional to n.
    A bound of Theta(n) is compared against the given categories.
    f(n) = Theta(n) is omega(log n), which corresponds to category 9.
    """
    # Case 1: max degree <= sqrt(log n). Worst case is a path P_n, steps = n-1.
    # n-1 is omega(log n), so category 9.
    f1_category = 9

    # Case 2: max degree <= log n. Worst case is a path P_n, steps = n-1.
    # n-1 is omega(log n), so category 9.
    f2_category = 9

    # Case 3: any forest. Worst case is a path P_n, steps = n-1.
    # n-1 is omega(log n), so category 9.
    f3_category = 9
    
    # The final answer is the three-digit number formed by these categories.
    final_answer = f"{f1_category}{f2_category}{f3_category}"
    print(final_answer)

solve()