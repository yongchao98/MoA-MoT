def solve():
    """
    This function analyzes the graph process and determines the categories for the bounds.
    The analysis is based on a construction that maximizes the number of steps, leading to a general bound
    on the number of steps S_max = Theta(min(Delta, sqrt(n))), where Delta is the maximum degree and
    n is the number of vertices.
    """

    # Case 1: Forest with max degree at most sqrt(log n)
    # S_max = Theta(min(sqrt(log n), sqrt(n))) = Theta(sqrt(log n))
    # We need to categorize sqrt(log n) = (log n)^0.5.
    # log(sqrt(log n)) = 0.5 * log(log n)
    # log(2^omega(sqrt(log log n))) = omega(sqrt(log log n))
    # Since 0.5*log(log n) is omega(sqrt(log log n)), it fits the lower bound of category 6.
    # The upper bound is O(log^0.9 n), which is also satisfied.
    f1_category = 6

    # Case 2: Forest with max degree at most log n
    # S_max = Theta(min(log n, sqrt(n))) = Theta(log n) for large n.
    # This directly corresponds to category 8.
    f2_category = 8

    # Case 3: Any forest
    # The bound is S_max = Theta(min(Delta, sqrt(n))).
    # The maximum number of steps is achieved when Delta is large, so S_max = Theta(sqrt(n)).
    # We need to categorize sqrt(n).
    # sqrt(n) is omega(log n), so it falls into category 9.
    f3_category = 9

    print(f"f1(n) for a forest of maximum degree at most sqrt(log n) falls into category: {f1_category}")
    print(f"f2(n) for a forest of maximum degree at most log n falls into category: {f2_category}")
    print(f"f3(n) for any forest falls into category: {f3_category}")

    final_number = f"{f1_category}{f2_category}{f3_category}"
    print(f"\nThe resulting three-digit number is: {final_number}")

solve()