def solve_graph_process():
    """
    Calculates the classification for the number of steps in a graph process for three cases.

    The analysis shows that the number of steps, T, is bounded by the maximum degree, Delta.
    Furthermore, a construction exists that achieves T = Theta(Delta) steps, as long as n >= O(Delta^2).

    Case 1: Maximum degree Delta <= sqrt(log n).
    The number of steps f1(n) is Theta(Delta) = Theta(sqrt(log n)).
    - sqrt(log n) is omega(2^sqrt(log log n)) and O((log n)^0.9).
    - This corresponds to category 6.
    """
    f1_category = 6

    """
    Case 2: Maximum degree Delta <= log n.
    The number of steps f2(n) is Theta(Delta) = Theta(log n).
    - We can construct a tree with Delta ~ log n lasting T ~ log n steps, as (log n)^2 <= n for large n.
    - Theta(log n) corresponds to category 8.
    """
    f2_category = 8

    """
    Case 3: Any forest.
    The number of steps T is maximized by a construction yielding T = Theta(sqrt(n)).
    - This requires Delta ~ sqrt(n) and n_vertices ~ Delta^2 ~ (sqrt(n))^2 ~ n.
    - f3(n) = Theta(sqrt(n)) is omega(log n).
    - This corresponds to category 9.
    """
    f3_category = 9

    # The final answer is the concatenation of the three category digits.
    final_answer_string = f"{f1_category}{f2_category}{f3_category}"
    
    print("The three digits for the bounds are derived as follows:")
    print(f"1) For f1(n) with max degree at most sqrt(log n), the bound is Theta(sqrt(log n)), which falls into category {f1_category}.")
    print(f"2) For f2(n) with max degree at most log n, the bound is Theta(log n), which falls into category {f2_category}.")
    print(f"3) For f3(n) on any forest, the bound is Theta(sqrt(n)), which falls into category {f3_category}.")
    
    print("\nThe resulting three-digit number is:")
    print(final_answer_string)


solve_graph_process()