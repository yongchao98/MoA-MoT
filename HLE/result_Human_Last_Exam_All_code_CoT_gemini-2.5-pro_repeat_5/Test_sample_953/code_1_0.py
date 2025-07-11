#
# This script determines the complexity classes for a variant of Luby's algorithm
# for finding a Maximal Independent Set on three different types of graphs.
#

def solve():
    """
    Analyzes the algorithm's complexity and determines the three-digit code.
    """

    # The algorithm is a randomized parallel MIS algorithm (Luby's algorithm).
    # We need to find the number of steps f(n) for it to finish with high probability.

    # --- Case 1: Cycle of length n ---
    # A cycle graph is 2-regular, so its maximum degree Delta is 2.
    # For graphs with specific structures like cycles (which are "tree-decomposable"),
    # the runtime of Luby's algorithm is known to be O(log Delta).
    # For a cycle, f1(n) = O(log 2) = O(1).
    # According to the problem's categories, O(1) corresponds to 1.
    d1 = 1

    # --- Case 2: Any tree on n vertices of degree at most 100 ---
    # Trees also have a structure for which the runtime is O(log Delta).
    # The maximum degree is bounded by Delta <= 100.
    # Therefore, f2(n) = O(log 100), which is a constant, O(1).
    # This also corresponds to category 1.
    d2 = 1

    # --- Case 3: Any graph on n vertices of degree at most 100 ---
    # For general graphs, even with a bounded degree, the runtime depends on n.
    # The standard analysis shows f3(n) = O(log n).
    # This bound is tight, as shown by expander graphs, which require Omega(log n) steps.
    # Thus, f3(n) = Theta(log n).
    # The category for functions that are Omega(log n) is 9.
    d3 = 9

    # The final answer is the concatenation of the digits d1, d2, and d3.
    # The prompt requests to output each number in the final equation.
    # We interpret this as showing the value for each digit.
    print(f"Analysis result for f1(n) [cycles]: Category {d1}")
    print(f"Analysis result for f2(n) [trees, deg<=100]: Category {d2}")
    print(f"Analysis result for f3(n) [general graphs, deg<=100]: Category {d3}")

    final_code = f"{d1}{d2}{d3}"
    print(f"\nThe final combined three-digit code is: {final_code}")

solve()