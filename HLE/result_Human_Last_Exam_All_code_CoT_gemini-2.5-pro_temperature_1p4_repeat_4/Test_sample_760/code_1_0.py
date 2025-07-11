import networkx as nx
import numpy as np

def check_graph_is_fixed_point(G):
    """
    Checks if a graph G satisfies the condition T(G) = G.
    A graph G is a fixed point if its adjacency matrix A is identical
    to the adjacency matrix of T(G).
    """
    num_vertices = G.number_of_nodes()
    
    # Graphs with 0 vertices are not considered. Graphs with 1 vertex are fixed points.
    if num_vertices == 0:
        return False
    if num_vertices == 1:
        return True

    # 1. Get the adjacency matrix A of the graph G.
    try:
        # Using a fixed nodelist ensures consistent matrix ordering.
        A = nx.to_numpy_array(G, nodelist=sorted(G.nodes()), dtype=np.int64)
    except Exception:
        # Fallback for older networkx versions
        A = nx.to_numpy_array(G, dtype=np.int64)

    # 2. Compute A^2. The entry (i, j) of A^2 gives the number of 
    #    length-2 paths between vertex i and vertex j.
    A_squared = np.dot(A, A)

    # 3. Construct the adjacency matrix of the transformed graph T(G).
    #    An edge (i,j) exists in T(G) if A_squared[i,j] is 1 or 2.
    A_transformed = np.zeros_like(A)
    condition_mask = (A_squared >= 1) & (A_squared <= 2)
    A_transformed[condition_mask] = 1

    # Ensure no self-loops in the transformed graph by setting the diagonal to 0.
    np.fill_diagonal(A_transformed, 0)
    
    # 4. Check if the original and transformed adjacency matrices are identical.
    return np.array_equal(A, A_transformed)

def find_all_fixed_point_graphs():
    """
    Finds all non-isomorphic, connected graphs with at most 7 vertices
    that are fixed points of the transformation T.
    """
    # We will store the count for each number of vertices from 1 to 7.
    counts_per_vertex_num = {n: 0 for n in range(1, 8)}

    # networkx.graph_atlas_g() provides all non-isomorphic graphs up to 7 vertices.
    all_graphs_in_atlas = nx.graph_atlas_g()

    for G in all_graphs_in_atlas:
        # The atlas list contains None for indices without graphs.
        if G is None:
            continue
            
        # We only consider connected graphs.
        if nx.is_connected(G):
            if check_graph_is_fixed_point(G):
                num_nodes = G.number_of_nodes()
                if num_nodes in counts_per_vertex_num:
                    counts_per_vertex_num[num_nodes] += 1
    
    # Print the results in a structured way.
    print("Number of qualifying graphs found for each vertex count:")
    total_count = 0
    equation_parts = []
    for n in range(1, 8):
        count_n = counts_per_vertex_num.get(n, 0)
        print(f"  - For {n} vertices: {count_n} graph(s)")
        total_count += count_n
        equation_parts.append(str(count_n))

    equation_str = " + ".join(equation_parts)
    print("\nThe total number of such graphs is the sum of the counts for each vertex size.")
    print(f"Total count = {equation_str} = {total_count}")

# Execute the main function to solve the problem.
find_all_fixed_point_graphs()