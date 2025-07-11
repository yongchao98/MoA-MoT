import networkx as nx
import numpy as np
from collections import defaultdict

def solve():
    """
    This script finds the number of non-isomorphic, connected graphs G with at most 7 vertices
    that satisfy the condition T(G) = G.

    The transformation T is defined as follows: an edge exists between distinct vertices x and y
    in T(G) if and only if there are exactly one or two paths of length 2 between x and y in G.

    The plan is to:
    1. Iterate through all non-isomorphic graphs up to 7 vertices using networkx.graph_atlas_g().
    2. For each graph, check if it's connected.
    3. If connected, check if it satisfies the T(G) = G condition using adjacency matrices.
       - A graph satisfies the condition if its adjacency matrix A is identical to the
         adjacency matrix A' derived from the transformation rule.
       - A'[i, j] = 1 if 1 <= (A^2)[i, j] <= 2, and 0 otherwise (for i != j).
    4. Count the graphs that satisfy the condition and report the total.
    """

    def check_graph_is_fixed_point(G):
        """
        Checks if a graph G satisfies the condition T(G) = G.
        """
        n = G.number_of_nodes()
        if n == 0:
            return False

        # Get the adjacency matrix A for the graph G.
        # The node order is sorted to ensure consistency.
        A = nx.to_numpy_array(G, nodelist=sorted(G.nodes()))

        # A^2[i, j] gives the number of length-2 paths between vertices i and j.
        A_squared = A @ A

        # Construct the adjacency matrix A_prime for the transformed graph T(G).
        A_prime = np.zeros((n, n), dtype=int)
        for i in range(n):
            for j in range(n):
                if i == j:
                    continue
                # Apply the transformation rule.
                if 1 <= A_squared[i, j] <= 2:
                    A_prime[i, j] = 1

        # The graph is a fixed point if its original adjacency matrix
        # is the same as the transformed one.
        return np.array_equal(A, A_prime)

    # Use a dictionary to store solution graphs, grouped by the number of vertices.
    solutions_by_n = defaultdict(list)
    
    # nx.graph_atlas_g() returns a list of all graphs with up to 7 vertices.
    all_graphs = nx.graph_atlas_g()

    # The graph's name in the atlas is G{i} where i is its index in the list.
    for i, G in enumerate(all_graphs):
        if G.number_of_nodes() == 0:
            continue
        
        # We only consider connected graphs.
        if nx.is_connected(G):
            # Check if the graph is a fixed point of the transformation T.
            if check_graph_is_fixed_point(G):
                graph_id = f"G{i}"
                solutions_by_n[G.number_of_nodes()].append(graph_id)

    total_count = 0
    count_per_n = []
    
    print("Found graphs satisfying the condition, broken down by number of vertices:")
    
    # Sort by n for a clean, ordered output.
    for n in sorted(solutions_by_n.keys()):
        graphs_found = solutions_by_n[n]
        count = len(graphs_found)
        total_count += count
        count_per_n.append(str(count))
        print(f" - For n={n}: {count} graph(s) found ({', '.join(graphs_found)})")
        
    print("-" * 30)
    if total_count > 0:
        # As requested, display the final count as a sum.
        equation = " + ".join(count_per_n)
        print(f"Final Equation: {equation} = {total_count}")
        print(f"The total number of such graphs is: {total_count}")
    else:
        print("No graphs found satisfying the condition.")

solve()