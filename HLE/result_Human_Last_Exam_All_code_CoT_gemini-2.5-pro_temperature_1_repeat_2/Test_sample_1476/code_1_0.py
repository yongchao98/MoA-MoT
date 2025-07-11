import numpy as np

def analyze_graph_signals():
    """
    This function demonstrates the logic from the problem description on an example graph.
    It shows that the given premises force the total variation of the vertex signal to be zero.
    """
    # Step 1: Define a graph with a cycle (a square) and a vertex signal x0.
    # To satisfy the problem's conditions, our analysis shows that the vertex signal x0
    # must be constant across all vertices in a connected graph.
    nodes = [0, 1, 2, 3]
    edges = [(0, 1), (1, 2), (2, 3), (3, 0)]
    constant_value = 10
    x0 = np.array([constant_value, constant_value, constant_value, constant_value])

    print("--- Analysis Setup ---")
    print(f"Graph: Vertices V={nodes}, Edges E={edges}")
    print(f"Vertex signal x0 = {x0}\n")

    # Step 2: Compute the edge signal x1 based on the definition x1_e = |x0_u - x0_v|.
    x1 = np.array([abs(x0[u] - x0[v]) for u, v in edges])

    print("--- Verifying Premises ---")
    print("Premise 1: 'No cycles having non-zero sum'.")
    # Since x1 is non-negative, this means the sum of x1 values on a cycle must be 0.
    # Our graph has one cycle involving all edges.
    cycle_sum = np.sum(x1)
    print(f"Sum of x1 values on the cycle {edges} = {cycle_sum}.")
    if cycle_sum == 0:
        print("Premise 1 is SATISFIED.\n")
    else:
        print("Premise 1 is NOT SATISFIED.\n")

    print("Premise 2: B1 * x1 = 0 (zero divergence).")
    # Define the vertex-edge incidence matrix B1 for the orientation 0->1, 1->2, 2->3, 3->0.
    B1 = np.array([
        [-1,  0,  0,  1],
        [ 1, -1,  0,  0],
        [ 0,  1, -1,  0],
        [ 0,  0,  1, -1]
    ])
    divergence = B1 @ x1
    print(f"The divergence vector B1 * x1 = {divergence.tolist()}.")
    if np.all(divergence == 0):
        print("Premise 2 is SATISFIED.\n")
    else:
        print("Premise 2 is NOT SATISFIED.\n")

    # Step 3: Draw the final conclusion about the Total Variation.
    # The conclusion is that if both premises hold, the Total Variation must be 0.
    total_variation = np.sum(x1)
    
    print("--- Final Inference ---")
    print("The premises together imply that the edge signal x1 must be the zero vector.")
    print("This means the Total Variation (TV) of the vertex signal x0 must be 0.")
    print("TV(x0) = sum over all edges {u,v} of |x0_u - x0_v|")

    # Display the equation for the total variation calculation.
    sum_str = " + ".join([f"|{x0[u]}-{x0[v]}|" for u,v in edges])
    val_str = " + ".join([f"{val}" for val in x1])

    print("\nFinal Calculation:")
    print(f"Total Variation = {sum_str}")
    print(f"                = {val_str}")
    print(f"                = {total_variation}")

# Run the analysis
analyze_graph_signals()