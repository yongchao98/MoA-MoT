import networkx as nx
import numpy as np

def solve():
    """
    Finds the number of non-isomorphic, connected graphs G with at most 7 vertices
    that satisfy the given transformation property T(G) = G.
    """
    
    def check_graph_T_equals_G(G):
        """
        Checks if a graph G satisfies the condition T(G) = G.

        T is a transformation where an edge (x,y) exists in T(G) if and only if
        there are exactly one or two paths of length 2 between x and y in G.
        """
        n = G.number_of_nodes()
        if n < 2:
            # A graph with 0 or 1 vertex has no pairs of distinct vertices.
            # T(G) will also have no edges, so T(G) = G is satisfied.
            return True

        # Use a sorted nodelist to get a canonical adjacency matrix representation.
        A = nx.to_numpy_array(G, nodelist=sorted(G.nodes()))

        # The number of length-2 paths between any two vertices is given by A^2.
        A_squared = np.dot(A, A)
        
        # We will construct the adjacency matrix for the transformed graph T(G).
        A_transformed = np.zeros_like(A)
        
        # An edge exists in T(G) if the corresponding entry in A_squared is 1 or 2.
        # This condition applies to pairs of distinct vertices.
        condition_mask = (A_squared >= 1) & (A_squared <= 2)
        A_transformed[condition_mask] = 1
        
        # Set diagonal to zero as we don't allow self-loops.
        np.fill_diagonal(A_transformed, 0)
        
        # The condition we are checking is T(G) = G, which means their 
        # adjacency matrices must be equal.
        return np.array_equal(A, A_transformed)

    counts_per_n = {}
    total_count = 0
    
    # Case n=1: The single-vertex graph. It is connected and satisfies the condition.
    counts_per_n[1] = 1
    
    # Case n=2: The only connected graph is an edge (K2). We check it.
    G2 = nx.complete_graph(2)
    if check_graph_T_equals_G(G2):
      counts_per_n[2] = 1
    else:
      counts_per_n[2] = 0

    # Cases n=3 to 7: We use the graph atlas from networkx.
    # The atlas contains all non-isomorphic graphs with 3 to 7 vertices.
    all_graphs_3_to_7 = nx.graph_atlas_g()
    
    for n in range(3, 8):
        count_for_n = 0
        # Filter the atlas for graphs with n vertices.
        graphs_for_n = [g for g in all_graphs_3_to_7 if g.number_of_nodes() == n]
        
        for G in graphs_for_n:
            # We are only interested in connected graphs.
            if nx.is_connected(G):
                if check_graph_T_equals_G(G):
                    count_for_n += 1
        
        counts_per_n[n] = count_for_n

    # Format the output as an equation showing the count for each n and the total.
    equation_parts = []
    for n in sorted(counts_per_n.keys()):
        total_count += counts_per_n[n]
        equation_parts.append(str(counts_per_n[n]))

    print("The number of such graphs for n=1, 2, 3, 4, 5, 6, 7 are, respectively:")
    print(" + ".join(equation_parts) + f" = {total_count}")
    print("\nThis means:")
    print(f"- For n=1, there is {counts_per_n[1]} graph.")
    print(f"- For n=2, there are {counts_per_n[2]} graphs.")
    print(f"- For n=3, there is {counts_per_n[3]} graph (the complete graph K3).")
    print(f"- For n=4, there is {counts_per_n[4]} graph (the complete graph K4).")
    print(f"- For n=5, there are {counts_per_n[5]} graphs.")
    print(f"- For n=6, there are {counts_per_n[6]} graphs.")
    print(f"- For n=7, there are {counts_per_n[7]} graphs.")
    
solve()
<<<3>>>