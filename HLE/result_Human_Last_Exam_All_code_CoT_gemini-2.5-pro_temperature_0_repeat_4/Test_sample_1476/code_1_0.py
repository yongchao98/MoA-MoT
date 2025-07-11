import networkx as nx
import numpy as np

def solve_and_explain():
    """
    This function demonstrates the logical deduction for the problem.
    It sets up a sample graph and signals that satisfy the premises,
    and then evaluates the given answer choices.
    """
    print("Step 1: Set up a graph G with a cycle and define signals.")
    # We use a square graph, which has a cycle.
    G = nx.cycle_graph(4)
    V = list(G.nodes())
    # We use an oriented edge list for the incidence matrix.
    E_oriented = [(0, 1), (1, 2), (2, 3), (3, 0)]
    E_unoriented = list(G.edges())

    print(f"Graph G has nodes V = {V} and edges E = {E_unoriented}")
    print(f"The graph has a cycle: {nx.cycle_basis(G)[0]}")

    # To satisfy the conditions, our analysis shows x1 must be 0.
    # This occurs if the vertex signal x0 is constant over the connected component.
    x0 = np.array([10.0, 10.0, 10.0, 10.0])
    print(f"\nLet's define a vertex signal x0 = {x0}")

    # Calculate x1 based on Condition 3.
    x1 = np.array([abs(x0[v] - x0[u]) for u, v in E_unoriented])
    print(f"This gives an edge signal x1 = {x1}")
    print("-" * 30)

    print("Step 2: Verify that the premises hold for this configuration.")
    # Premise 1: No cycles with non-zero sum.
    # Since x1 is all zeros, the sum of its values over any cycle is 0.
    print("Premise 1 (no cycles with non-zero sum): Satisfied, as all x1 values are 0.")

    # Premise 2: B1 * x1 = 0.
    B1 = nx.incidence_matrix(G, nodelist=V, edgelist=E_oriented, oriented=True).toarray()
    divergence = B1 @ x1
    print(f"Premise 2 (B1 @ x1 = 0): Satisfied. The divergence is {divergence}.")
    print("-" * 30)

    print("Step 3: Evaluate the answer choices based on our deduction that x1 must be 0.")
    # A. x1 is an eigenvector of L1.
    print("A. x1 is an eigenvector of L1: False. The zero vector is by definition not an eigenvector.")

    # B. The graph G has no cycles.
    print(f"B. The graph G has no cycles: False. Our example graph G has a cycle.")

    # C. x1 is in ker(B1).
    print(f"C. x1 is in ker(B1): True. This is equivalent to Premise 2 (B1 @ x1 = 0), which holds.")

    # D. G has total variation = 0.
    total_variation = np.sum(x1)
    x1_list = list(x1)
    sum_str = " + ".join(f"{val:.1f}" for val in x1_list)
    print(f"D. G has total variation = 0: True.")
    print(f"   Total Variation = sum(|x0_u - x0_v|) = sum(x1)")
    print(f"   Equation: {sum_str} = {total_variation:.1f}")
    print("-" * 30)

    print("Step 4: Final Conclusion.")
    print("The premises force the edge signal x1 to be the zero vector.")
    print("This implies the total variation of the graph signal is 0 (Choice D).")
    print("While Choice C is also true, it is a weaker statement (a direct premise).")
    print("Choice D is the strongest conclusion derived from all premises combined.")

solve_and_explain()