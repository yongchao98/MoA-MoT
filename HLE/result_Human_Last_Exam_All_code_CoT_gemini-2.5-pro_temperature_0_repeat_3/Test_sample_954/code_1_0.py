def solve():
    """
    This function encapsulates the reasoning for determining the bounds on the simulation process.
    """

    # Step 1 & 2: Analysis and Upper Bound
    # The number of steps T is bounded by the maximum degree Delta of the forest.
    # In each round, any vertex 'u' with the current maximum degree 'd_max' loses 1 life point.
    # Loss(u) = sum(1/max(d(u), d(v)) for v in neighbors(u))
    # Since d(u) = d_max, max(d(u), d(v)) is always d_max.
    # Loss(u) = sum(1/d_max) = d_max * (1/d_max) = 1.
    # Since vertices start with 1 life, all max-degree vertices die in the first round.
    # The max degree of the remaining graph is strictly smaller.
    # This implies the process lasts at most Delta rounds. So, T = O(Delta).

    # Step 3: Lower Bound
    # A construction of a tree with a path of vertices v_1, ..., v_Delta, where
    # deg(v_i) = i (achieved by adding leaves), shows that the process can take
    # Omega(Delta) steps. In this tree, only the highest-degree vertices are removed
    # in each step, peeling the tree one layer at a time.
    # Therefore, the number of steps T = Theta(Delta).

    # Step 4: Apply to each case
    
    # Case 1: Forest with max degree Delta <= sqrt(log n)
    # The maximum number of steps is T = Theta(sqrt(log n)) = Theta((log n)^0.5).
    # We check the given categories for f(n) = (log n)^0.5.
    # Category 6 is: f(n) = 2^omega(sqrt(log log n)) but f(n) = O((log n)^0.9).
    # - Check lower bound: log(f(n)) = 0.5 * log(log n).
    #   log(2^sqrt(log log n)) = sqrt(log log n).
    #   Since 0.5 * log(log n) is omega(sqrt(log log n)), the lower bound holds.
    # - Check upper bound: (log n)^0.5 is O((log n)^0.9).
    # So, the first digit is 6.
    f1_category = 6

    # Case 2: Forest with max degree Delta <= log n
    # The maximum number of steps is T = Theta(log n).
    # This directly matches category 8: f(n) = Theta(log n).
    # So, the second digit is 8.
    f2_category = 8

    # Case 3: Any forest
    # The max degree Delta can be up to n-1 (for a star graph).
    # The maximum number of steps is T = Theta(Delta) = Theta(n).
    # We check the categories for f(n) = Theta(n).
    # Theta(n) is omega(log n). This fits category 9.
    # So, the third digit is 9.
    f3_category = 9

    # The final three-digit number is formed by these categories.
    final_number_str = f"{f1_category}{f2_category}{f3_category}"
    
    print(f"Analysis result:")
    print(f"1) For a forest of maximum degree at most sqrt(log n), the bound f1(n) corresponds to category: {f1_category}")
    print(f"2) For a forest of maximum degree at most log n, the bound f2(n) corresponds to category: {f2_category}")
    print(f"3) For any forest, the bound f3(n) corresponds to category: {f3_category}")
    print(f"\nThe resulting three-digit number is: {final_number_str}")

solve()