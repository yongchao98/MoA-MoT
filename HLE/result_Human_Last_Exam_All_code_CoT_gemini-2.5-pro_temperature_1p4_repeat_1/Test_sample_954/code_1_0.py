def solve():
    """
    This function encapsulates the reasoning for determining the bounds on the number of steps
    for the described graph erosion process.
    """

    # Case 3: On any forest.
    # The literature on this process shows that for some tree structures, the number of steps
    # can be Theta(log n), and it's bounded by O(log n) for all trees.
    # Therefore, the maximum number of steps is Theta(log n).
    # This corresponds to option 8.
    f3 = 8

    # Case 2: On any forest of maximum degree at most log n.
    # The constructions that yield a Theta(log n) runtime require a maximum degree
    # of Theta(log n). Since this is permitted by the constraint Delta <= log n,
    # the maximum number of steps is also Theta(log n).
    # This corresponds to option 8.
    f2 = 8

    # Case 1: On any forest of maximum degree at most sqrt(log n).
    # The constructions for Theta(log n) are not possible under this stricter degree constraint.
    # The tightest known bound for this case is T = O(log n / log Delta).
    # This bound is maximized for the largest possible Delta in this class, i.e., Delta = sqrt(log n).
    # The complexity is O(log n / log(sqrt(log n))) = O(log n / log(log n)).
    # Let's compare this to the provided options:
    # - Is it omega(log^{0.9} n)? Yes, because log^{0.1} n grows faster than log(log n).
    # - Is it o(log n)? Yes, because (log n / log log n) / log n = 1 / log log n, which tends to 0.
    # This complexity class matches option 7.
    f1 = 7
    
    # The problem asks for a three-digit number where the i-th digit corresponds to f_i(n).
    final_number_string = f"{f1}{f2}{f3}"
    print(final_number_string)

solve()