import math

def solve():
    """
    Analyzes the runtime of a randomized MIS algorithm and provides the complexity class for three graph families.
    """

    # Analysis of the Algorithm
    # The algorithm is a variant of Luby's algorithm for Maximal Independent Set (MIS).
    # In each step, every vertex v gets a random value r(v) in [0,1]. If r(v) is greater than r(w) for all neighbors w, v is added to the MIS. Then v and all its neighbors are removed.
    # We want to find the number of steps f(n) for the algorithm to terminate w.h.p. (with probability 1 - 1/n).

    # Step 1: General Upper Bound (Case 3, and thus Cases 1 & 2)
    # Let G be a graph with maximum degree at most Δ = 100.
    # Consider an arbitrary edge e = (u, v) in the graph at some step.
    # The vertex u is selected if its random value is the highest among itself and its neighbors. Let d(u) be its current degree. The probability of this is P(u is selected) = 1 / (d(u) + 1).
    # Similarly, P(v is selected) = 1 / (d(v) + 1).
    # The events "u is selected" and "v is selected" are disjoint, because the former implies r(u) > r(v) and the latter implies r(v) > r(u).
    # If either u or v is selected, the edge (u, v) is removed from the graph.
    # So, P(edge (u,v) is removed) >= P(u is selected or v is selected) = P(u is selected) + P(v is selected).
    # Since d(u) <= 100 and d(v) <= 100, we have P(edge (u,v) is removed) >= 1/(100+1) + 1/(100+1) = 2/101.
    # This means every edge has at least a constant probability of being removed in each step.
    # In expectation, a constant fraction (at least 2/101) of the edges are removed in each step.
    # Concentration bounds show that the number of edges removed is close to its expectation w.h.p.
    # Thus, the number of edges decreases geometrically. To reduce the number of edges from O(n) to zero, it takes O(log n) steps.
    # This O(log n) bound applies to any graph with degree at most 100.
    # Therefore, f_3(n) = O(log n). Since cycles and trees are specific types of these graphs, f_1(n) = O(log n) and f_2(n) = O(log n) as well.

    # Step 2: Lower Bounds for each case
    # To find the best possible function (tight bound), we check for a lower bound.
    # We need to find if there's a graph in each class that forces Ω(log n) steps.

    # Case 1: A cycle of length n.
    # A cycle C_n has a maximum degree of 2. It's a member of the general class of graphs with Δ <= 100.
    # It is a standard result in the analysis of parallel MIS algorithms that on a path or a cycle, this algorithm takes Ω(log n) steps w.h.p.
    # Intuitively, in each step, only a fraction of vertices are removed, leaving a collection of smaller paths/cycles. This process must be repeated logarithmically many times.
    # Since f_1(n) = O(log n) and f_1(n) = Ω(log n), we have f_1(n) = Θ(log n).
    # According to the provided categories, this corresponds to 9: f(n) = Ω(log n).
    d1 = 9

    # Case 2: Any tree on n vertices of degree at most 100.
    # The upper bound is f_2(n) = O(log n).
    # For the lower bound, we must consider the worst-case tree. A simple path P_n is a tree with a maximum degree of 2, so it belongs to this class.
    # As mentioned, the algorithm takes Ω(log n) steps on a path.
    # Thus, the worst-case performance for this class of trees is Ω(log n).
    # Since f_2(n) = O(log n) and f_2(n) = Ω(log n), we have f_2(n) = Θ(log n).
    # This also corresponds to category 9.
    d2 = 9

    # Case 3: Any graph on n vertices of degree at most 100.
    # The upper bound is f_3(n) = O(log n).
    # For the lower bound, we note that this class includes cycles and paths.
    # Since the algorithm takes Ω(log n) on these specific graphs, the worst-case performance for this general class must also be Ω(log n).
    # Since f_3(n) = O(log n) and f_3(n) = Ω(log n), we have f_3(n) = Θ(log n).
    # This also corresponds to category 9.
    d3 = 9

    # Step 3: Final Answer
    # We combine the digits for each case.
    final_answer = f"{d1}{d2}{d3}"

    print("Analysis of the algorithm's runtime:")
    print("1. For a cycle of length n, f_1(n) = Theta(log n).")
    print("2. For any tree with max degree <= 100, the worst-case runtime is on a path, so f_2(n) = Theta(log n).")
    print("3. For any graph with max degree <= 100, the worst-case runtime is on a path or cycle, so f_3(n) = Theta(log n).")
    print("")
    print("Mapping to categories:")
    print(f"f_1(n) = Theta(log n) falls into category {d1}: f(n) = Omega(log n).")
    print(f"f_2(n) = Theta(log n) falls into category {d2}: f(n) = Omega(log n).")
    print(f"f_3(n) = Theta(log n) falls into category {d3}: f(n) = Omega(log n).")
    print("")
    print(f"The resulting three-digit number d1 d2 d3 is: {final_answer}")
    print("<<<999>>>")

solve()