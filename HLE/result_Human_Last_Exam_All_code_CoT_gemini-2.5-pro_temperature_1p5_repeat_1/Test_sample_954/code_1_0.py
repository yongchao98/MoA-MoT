def solve():
    """
    This function provides the final result based on the analysis of the graph process.
    The analysis concludes that the process terminates in a constant number of steps
    for all types of forests described.
    """
    # Case 1: f_1(n) on any forest of maximum degree at most sqrt(log n).
    # The number of steps is O(1). This falls into category 1.
    f1_category = 1

    # Case 2: f_2(n) on any forest of maximum degree at most log n.
    # The number of steps is O(1). This falls into category 1.
    f2_category = 1

    # Case 3: f_3(n) on any forest.
    # The number of steps is O(1). This falls into category 1.
    f3_category = 1

    # The final result is the concatenation of the category digits.
    final_answer = f"{f1_category}{f2_category}{f3_category}"
    print(final_answer)

solve()