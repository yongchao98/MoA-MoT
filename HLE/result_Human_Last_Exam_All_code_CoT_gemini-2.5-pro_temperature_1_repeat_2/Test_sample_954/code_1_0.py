def solve():
    """
    This function encapsulates the reasoning to determine the three-digit number.
    The reasoning is based on analyzing the "Bullying Process" on forests.

    Case 1: Maximum degree at most sqrt(log n).
    The number of steps T is bounded by the maximum degree Delta.
    T <= Delta <= sqrt(log n).
    The function f(n) = sqrt(log n) is omega(2^sqrt(log log n)) and O(log^0.9 n).
    This corresponds to category 6.

    Case 2: Maximum degree at most log n.
    It's possible to construct a tree where the number of steps T is Theta(log n)
    while the maximum degree Delta is also Theta(log n). This falls within the constraint.
    f(n) = Theta(log n) corresponds to category 8.

    Case 3: Any forest.
    The maximum number of steps for any forest is Theta(log n).
    f(n) = Theta(log n) corresponds to category 8.

    The resulting three-digit number is 688.
    """
    f1_category = 6
    f2_category = 8
    f3_category = 8
    
    final_answer = f1_category * 100 + f2_category * 10 + f3_category
    print(final_answer)

solve()