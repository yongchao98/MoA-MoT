import math

def solve():
    """
    Analyzes the runtime of a randomized MIS algorithm on different graph classes
    and determines the corresponding complexity category.
    """

    # Define the categories for f(n)
    categories = {
        1: "O(1)",
        2: "omega(1) but o(log*n)",
        3: "Theta(log*n)",
        4: "omega(log*n) but o(log log n)",
        5: "Theta(log log n)",
        6: "omega(log log n) but O(log^0.1(n))",
        7: "omega(log^0.1(n)) but O(log^0.9(n))",
        8: "omega(log^0.9(n)) but o(log n)",
        9: "Omega(log n)",
    }

    # --- Analysis for f_1(n) on a cycle ---
    # In the first step, the algorithm selects vertices that are local maxima based on their random values.
    # This breaks the cycle into a collection of disjoint paths.
    # It's a known result that for a random permutation on a cycle, the longest segment that
    # doesn't contain a local maximum has length O(log n) with high probability.
    # The runtime of the algorithm is then dominated by the time to find an MIS on the longest
    # of these paths. The runtime for a path of length k is Theta(log k).
    # Therefore, the total time is 1 + Theta(log(O(log n))) = Theta(log log n).
    d1 = 5
    print("Analysis for f_1(n) (Cycles):")
    print(f"The algorithm takes Theta(log(log n)) steps. This corresponds to category {d1}: {categories[d1]}.\n")

    # --- Analysis for f_2(n) on a tree of degree at most 100 ---
    # For any graph with maximum degree Delta, the algorithm terminates in O(log n) steps
    # with high probability. This is because in each step, any given vertex is removed with
    # a probability bounded below by a constant depending on Delta. This gives an O(log n) upper bound.
    # A path graph is a tree with a maximum degree of 2. For a path, the algorithm
    # requires Omega(log n) steps in the worst case.
    # Since the class of trees includes path graphs, the lower bound is Omega(log n).
    # Combining the upper and lower bounds, the complexity is Theta(log n).
    d2 = 9
    print("Analysis for f_2(n) (Trees with degree <= 100):")
    print(f"The algorithm takes Theta(log n) steps. This corresponds to category {d2}: {categories[d2]}.\n")

    # --- Analysis for f_3(n) on any graph of degree at most 100 ---
    # The O(log n) upper bound for any graph with bounded maximum degree applies here as well.
    # The class of general graphs with degree at most 100 also includes path graphs,
    # so the Omega(log n) lower bound also applies.
    # Therefore, the complexity is Theta(log n).
    d3 = 9
    print("Analysis for f_3(n) (Graphs with degree <= 100):")
    print(f"The algorithm takes Theta(log n) steps. This corresponds to category {d3}: {categories[d3]}.\n")

    # --- Final Answer ---
    final_answer_str = f"{d1}{d2}{d3}"
    print("The three digits are:")
    print(f"d1 = {d1}")
    print(f"d2 = {d2}")
    print(f"d3 = {d3}")
    print(f"\nThe final encoded answer is: {final_answer_str}")

solve()