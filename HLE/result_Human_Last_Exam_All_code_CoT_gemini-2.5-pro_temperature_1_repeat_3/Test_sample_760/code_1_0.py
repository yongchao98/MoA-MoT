import networkx as nx
import numpy as np

def solve():
    """
    Finds the number of non-isomorphic, connected graphs G with at most 7 vertices
    such that T(G) = G, where T is the specified transformation.
    """
    count = 0
    found_graphs_details = []

    # The networkx graph atlas contains all graphs with up to 7 vertices.
    # We iterate through them.
    for G in nx.graph_atlas_g():
        n = G.number_of_nodes()

        # We only consider connected graphs with at least one vertex.
        if n > 0 and nx.is_connected(G):
            # Get the adjacency matrix A of the graph G.
            # We sort the nodes to ensure a consistent ordering.
            nodelist = sorted(G.nodes())
            A = nx.to_numpy_array(G, nodelist=nodelist, dtype=int)

            # A^2[i,j] gives the number of length-2 paths between i and j.
            A_squared = np.dot(A, A)

            # Construct the adjacency matrix A_T for the transformed graph T(G).
            A_T = np.zeros_like(A, dtype=int)
            is_fixed_point = True
            
            # For each pair of distinct vertices (i, j), check the transformation condition.
            for i in range(n):
                for j in range(i + 1, n):
                    num_paths = A_squared[i, j]
                    
                    # An edge exists in T(G) iff there are 1 or 2 paths of length 2.
                    edge_in_T = 1 if (num_paths == 1 or num_paths == 2) else 0
                    
                    # If T(G) is different from G, it's not a fixed point.
                    if A[i, j] != edge_in_T:
                        is_fixed_point = False
                        break
                if not is_fixed_point:
                    break
            
            # If the loop completes and the graph is a fixed point, count it.
            if is_fixed_point:
                count += 1
                num_edges = G.number_of_edges()
                detail = f"Found a graph with {n} vertices and {num_edges} edges."
                found_graphs_details.append(detail)
                print(detail)

    print("\nSummary of found graphs:")
    # We found one graph for each of the cases detailed above.
    # The final count is the sum of these individual finds.
    equation_str = " + ".join(["1"] * count) + f" = {count}"
    print(f"The final count is derived from the sum: {equation_str}")
    
    print(f"\nTotal number of graphs satisfying the condition: {count}")

solve()
<<<5>>>